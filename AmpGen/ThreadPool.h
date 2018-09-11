#ifndef AMPGEN_THREADPOOL_H
#define AMPGEN_THREADPOOL_H

/* @class ThreadPool ThreadPool.h AmpGen/ThreadPool.h
 * Thread pool implementation taken from https://github.com/progschj/ThreadPool
 * Modified to allow explicit clearing of queues.
 * A single static thread pool exists that can be used as a sceduler.
 */

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace AmpGen
{
  class ThreadPool
  {
    public:
      static size_t nThreads;

      ~ThreadPool();

      template <class F, class... Args>
        static auto schedule( F&& f, Args&&... args ) -> 
        std::future<typename std::result_of<F( Args... )>::type>
        {
          return getMe()->enqueue( std::forward<F>( f ), std::forward<Args>( args )... );
        }


      template<class F, class... Args>
        auto enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(Args...)>::type>;
      ThreadPool(const size_t& nt);
    private:
      static ThreadPool* gThreadPool;
      static ThreadPool* getMe()
      {
        if ( !gThreadPool ) 
          gThreadPool = new ThreadPool(ThreadPool::nThreads);
        return gThreadPool;
      }

      std::vector<std::thread>          m_workers;
      std::queue<std::function<void()>> m_tasks;
      std::mutex                        m_queue_mutex;
      std::condition_variable           m_condition;
      bool                              m_stop;
  };

  template<class F, class... Args>
    auto ThreadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::result_of<F(Args...)>::type>
    {
      using return_type = typename std::result_of<F(Args...)>::type;

      auto task = std::make_shared< std::packaged_task<return_type()> >(
          std::bind(std::forward<F>(f), std::forward<Args>(args)...)
          );

      std::future<return_type> res = task->get_future();
      {
        std::unique_lock<std::mutex> lock(m_queue_mutex);

        if(m_stop)
          throw std::runtime_error("enqueue on stopped ThreadPool");
        m_tasks.emplace([task](){ (*task)(); });
      }
      m_condition.notify_one();
      return res;
    }
} //namespace AmpGen;
#endif
