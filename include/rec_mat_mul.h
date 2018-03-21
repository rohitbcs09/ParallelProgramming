#ifndef _REC_MAT_MUL_H
#define _REC_MAT_MUL_H

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <deque>
#include <stack>
#include <stdint.h>

using namespace std;
typedef std::vector<std::vector<uint64_t> > Matrix;

typedef void (*tasklet)(Matrix& , Matrix& , Matrix& , int, int, int, int, int, 
                        int, int, int);

void print(Matrix& a1, Matrix& a2, Matrix& a3, int x_row, int x_col, int y_row, 
           int y_col, int z_row, int z_col, int n, int id) {
    for(int i = 0; i<4; ++i) {
        for(int j = 0; j<4; ++j) {
            std::cout << a1[i][j]<<" ";
        }
        std::cout << "\n";
    }
}

typedef struct tasklet_desc {
  Matrix* X_;
  Matrix* Y_;
  Matrix* Z_;
  int x_row;
  int x_col;
  int y_row;
  int y_col;
  int z_row;
  int z_col;
  int n_;
} tasklet_desc;

typedef struct work {
    tasklet t_;
    tasklet_desc td_;
} Work;


typedef struct sync {
    Matrix* X_;
    Matrix* Y_;
    Matrix* Z_11;
    Matrix* Z_12;
    Matrix* Z_21;
    Matrix* Z_22;
    int quad;
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
        int count = 0;
        while(1) {
            if(!d_.empty()) {
                std::cout << "Inside Run : Popping Job" << "\n";
                Work w = this->pop_back();
                (*(w.t_))(*(w.td_.X_), *(w.td_.Y_), *(w.td_.Z_), w.td_.x_row,
                          w.td_.x_col, w.td_.y_row, w.td_.y_col, w.td_.z_row, 
                          w.td_.z_col, w.td_.n_, id_);
            }
            else {
                ++count;
                if(count == 20) {
                    std::cout << "Inside Run : Steal Failure after 20 attempts" << "\n";
                    return;
                }
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

#endif
