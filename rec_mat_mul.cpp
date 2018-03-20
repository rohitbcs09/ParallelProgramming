#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <deque>
#include <stack>
#include <stdint.h>

using namespace std;

typedef std::vector<std::vector<uint64_t> >
        (*tasklet)(std::vector<std::vector<uint64_t> > ,
                   std::vector<std::vector<uint64_t> > , int);

std::vector<std::vector<uint64_t> > 
    print(std::vector<std::vector<uint64_t> > a1,
          std::vector<std::vector<uint64_t> > a2, int n) {

    std::vector<std::vector<uint64_t> > res;
    for(int i = 0; i<4; ++i) {
        for(int j = 0; j<4; ++j) {
            std::cout << a1[i][j]<<" ";
        }
        std::cout << "\n";
    }
    return res;
}

typedef struct tasklet_desc {
  std::vector<std::vector<uint64_t> >* X_;
  std::vector<std::vector<uint64_t> >* Y_;
  std::vector<std::vector<uint64_t> >* Z_;
  int n_;
} tasklet_desc;

typedef struct work {
    tasklet t_;
    tasklet_desc td_;
} Work;

typedef struct sync {
    std::vector<std::vector<uint64_t> >* X_;
    std::vector<std::vector<uint64_t> >* Y_;
    std::vector<std::vector<uint64_t> >* Z_;
    int type_;
    int val_;
} sync;


class Task { 
public:
    Task(int id) {
        isDone = false;
        start = false;
        id_ = id;
    }
    
    ~Task() {
        d_.clear();
    }

    void push_front(Work w) {
        std::unique_lock<std::mutex> lock(m_);
        d_.push_front(w);
    }

    void push_back(Work w) {
        std::unique_lock<std::mutex> lock(m_);
        d_.push_back(w);
    }

    Work pop_back() {
        std::unique_lock<std::mutex> lock(m_);
        Work ret = d_.back();
        d_.pop_back();
        return ret;
    }

    Work pop_front() {
        std::unique_lock<std::mutex> lock(m_);
        Work ret = d_.front();
        d_.pop_front();
        return ret;
    }

    void run() {
        while(1) {
            if(!d_.empty()) {
                Work w = this->pop_back();
                (*(w.t_))(*(w.td_.X_), *(w.td_.Y_), w.td_.n_);
                return;
            }
            else {
                // STEAL
            }
        }
        return;
    }

private:
    std::deque<Work> d_;
    std::stack<sync> s_;
    std::mutex m_;
    bool start;
    bool isDone;
    int id_;
};

int main() {

    std::vector<Task *> pool;

    std::vector< std::vector<uint64_t> > A(4, std::vector<uint64_t>(4, 1)); 
    std::vector< std::vector<uint64_t> > B(4, std::vector<uint64_t>(4, 2)); 
    std::vector< std::vector<uint64_t> > Z(4, std::vector<uint64_t>(4, 0));

    Work work;
    work.t_ = &print;
    work.td_.X_ = &A;
    work.td_.Y_ = &B;
    work.td_.Z_ = &Z;
    work.td_.n_ = 4;

    unsigned int num_threads = std::thread::hardware_concurrency();
    for(int i = 0; i<num_threads; ++i) {
        pool.push_back(new Task(i));
    }
    
    // cilk spawn each thread in the pool
    for(int i = 0; i<num_threads; ++i) {
        pool[i]->push_back(work);
        pool[i]->run();
    }
    
    return 1;
}
