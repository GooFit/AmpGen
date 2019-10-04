#include "AmpGen/ThreadPool.h"

#include <iostream>

#include "AmpGen/MsgService.h"

using namespace AmpGen;

ThreadPool::ThreadPool(const size_t& nt)
{ 
  for ( size_t i = 0; i < nt; ++i )
    m_workers.emplace_back( [this] {
     for(;;){
        std::function<void()> task;
        {
          std::unique_lock<std::mutex> lock( this->m_queue_mutex );
          this->m_condition.wait( lock, [this] { return this->m_stop || !this->m_tasks.empty(); } );
          if( this->m_stop && this->m_tasks.empty()) return; 
          task = std::move( this->m_tasks.front() ); 
          this->m_tasks.pop();
        }
        task();  
      }
    } );
}

ThreadPool::~ThreadPool()
{
  waitForStoppedThreads();
}

void ThreadPool::waitForStoppedThreads()
{
  {
    std::unique_lock<std::mutex> lock(m_queue_mutex);
    m_stop = true;
  }
  m_condition.notify_all();
  for(std::thread &worker: m_workers) worker.join();
  m_stop = false;
}
