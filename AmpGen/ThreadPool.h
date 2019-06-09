#ifndef AMPGEN_THREADPOOL_H
#define AMPGEN_THREADPOOL_H

/* @class ThreadPool ThreadPool.h AmpGen/ThreadPool.h
 * Thread pool implementation taken from https://github.com/progschj/ThreadPool
 * Modified to allow explicit clearing of queues.
 *
 * Copyright (c) 2012 Jakob Progsch, Václav Zeman
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgment in the product documentation would be
 * appreciated but is not required.
 *
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 *
 * 3. This notice may not be removed or altered from any source distribution.
 */

#include <memory.h>
#include <stddef.h>
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
      ThreadPool(const size_t& nt);
      ~ThreadPool();
      template<class F, class... Args> auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>;
      void waitForStoppedThreads(); 

    private:
      std::vector<std::thread>          m_workers;
      std::queue<std::function<void()>> m_tasks;
      std::mutex                        m_queue_mutex;
      std::condition_variable           m_condition;
      bool                              m_stop={false};
  };

  template<class F, class... Args> auto ThreadPool::enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
  {
    using return_type = typename std::result_of<F(Args...)>::type;
    auto task = std::make_shared< std::packaged_task<return_type()> >( f, args... );
    std::future<return_type> res = task->get_future();
    {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      if(m_stop) throw std::runtime_error("enqueue on stopped ThreadPool");
      m_tasks.emplace([task](){ (*task)(); });
    }
    m_condition.notify_one();
    return res;
  }
} //namespace AmpGen;
#endif
