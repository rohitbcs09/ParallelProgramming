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

typedef struct sync {
    int type_;
    int val_;
    int x_row;
    int x_col;
    int y_row;
    int y_col;
    int z_row;
    int z_col;
    int n;
    struct sync* parent;

    sync(int type) {
        type_ = type;
        val_  = 4;
    }
} Sync;

void set_sync_state(Sync* c_sync, Sync *p_sync, int type, int x1, int x2, int y1, int y2, int z1, int z2, int size) {
    c_sync->type_ = type;
    c_sync->type_ = 4;
    c_sync->x_row = x1;
    c_sync->x_col = x2;
    c_sync->y_row = y1;
    c_sync->y_col = y2;
	c_sync->z_row = z1;
	c_sync->z_col = z2;
	c_sync->n = size;
	c_sync->parent = p_sync;
}

typedef void (*tasklet)(Matrix* , Matrix* , Matrix* , int, int, int, int, int, 
                        int, int, Sync*, int);

void print(Matrix* a1, Matrix* a2, Matrix* a3, int x_row, int x_col, int y_row, 
           int y_col, int z_row, int z_col, int n, Sync* s, int id) {
    for(int i = 0; i<4; ++i) {
        for(int j = 0; j<4; ++j) {
            std::cout << (*a1)[i][j]<<" ";
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
    Sync* sync_;
} Work;

class Task { 
public:
    Task(int id) {
        is_done = false;
        is_started = false;
        id_ = id;
    }
    
    ~Task() {
        d_.clear();
    }
    
    void set_is_started(bool val) {
        is_started = val;
    }

    void set_is_done(bool val) {
        is_done = val;
    }

    void push_front(Work* w) {
        std::unique_lock<std::mutex> lock(m_);
        d_.push_front(w);
    }

    void push_back(Work* w) {
        std::unique_lock<std::mutex> lock(m_);
        d_.push_back(w);
    }

    Work* pop_back() {
        std::unique_lock<std::mutex> lock(m_);
        Work *ret = d_.back();
        d_.pop_back();
        return ret;
    }

    Work* pop_front() {
        std::unique_lock<std::mutex> lock(m_);
        Work *ret = d_.front();
        d_.pop_front();
        return ret;
    }

    int get_deque_size() {
        return d_.size();
    }
    
    void push_sync(Sync* s) {
	    if(s == NULL) {
		    return;
		}
        std::unique_lock<std::mutex> lock(s_m_);
        s_.push(s);
    }

    Sync stack_top() {
        std::unique_lock<std::mutex> lock(s_m_);
        Sync *ret = s_.top();
        return *ret;
    }

    Sync* stack_top_ref() {
        std::unique_lock<std::mutex> lock(s_m_);
        Sync *ret = s_.top();
        return ret;
    }

    void dec_sync_ref(Sync *s) {
        std::unique_lock<std::mutex> lock(s_m_);
 	    if(s == NULL) {
		    return;
		}
        --s->val_; 
    }

    bool is_sync_ref_zero(Sync *s) {
 	    std::unique_lock<std::mutex> lock(s_m_);
        if(s == NULL) {
		    return false;
		}
        return (s->val_ == 0) ? true : false; 
    }
   
    int get_sync_ref_count(Sync *s) {
        std::unique_lock<std::mutex> lock(s_m_);
 	    if(s == NULL) {
		    return -1;
		}
        return s->val_;
    }

    void pop_sync() {
        std::unique_lock<std::mutex> lock(s_m_);
	    if(s_.size() == 0) {
		    return;
		}
        s_.pop();
    }
    
    bool is_stack_empty() {
        return s_.size() == 0 ? true : false;
    }

    void run() {
        int count = 0;
        while(1) {
            if(!d_.empty()) {
			    is_started = true;
                //std::cout << "Inside Run : Popping Job" << "\n";
                Work *w = this->pop_back();
                (*(w->t_))(w->td_.X_, w->td_.Y_, w->td_.Z_, w->td_.x_row,
                          w->td_.x_col, w->td_.y_row, w->td_.y_col, w->td_.z_row, 
                          w->td_.z_col, w->td_.n_, w->sync_, id_);
            }
            else {
                ++count;
                if(count == 500) {
                    std::cout << "Inside Run : Steal Failure after 20 attempts" << "\n";
                    return;
                }
                // STEAL
            }
        }
        return;
    }

private:
    std::deque<Work *> d_;
    std::stack<Sync *> s_;
    std::mutex m_;
    std::mutex s_m_;
    bool is_started;
    bool is_done;
    int id_;
};

#endif

