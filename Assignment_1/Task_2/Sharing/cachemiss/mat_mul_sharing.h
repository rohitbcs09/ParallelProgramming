#ifndef _MAT_MUL_STEAL_H
#define _MAT_MUL_STEAL_H

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <deque>
#include <stack>
#include <stdint.h>
#include <math.h> 

using namespace std;

extern uint64_t g_seed;

uint64_t fastrand() {
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
} 

typedef std::vector<std::vector<uint64_t> > Matrix;

typedef struct sync {
    int s_type;
    int s_val;
    int s_x_row;
    int s_x_col;
    int s_y_row;
    int s_y_col;
    int s_z_row;
    int s_z_col;
    int s_n;
    struct sync* s_parent;
} Sync;

Sync* create_sync_state(Sync *p_sync, int type, int x1, int x2, int y1, int y2, 
                        int z1, int z2, int size) {
    Sync *c_sync = new Sync;
    c_sync->s_type = type;
    c_sync->s_val = 4;
    c_sync->s_x_row = x1;
    c_sync->s_x_col = x2;
    c_sync->s_y_row = y1;
    c_sync->s_y_col = y2;
    c_sync->s_z_row = z1;
    c_sync->s_z_col = z2;
    c_sync->s_n = size;
    c_sync->s_parent = p_sync;
    return c_sync;
}

typedef void (*tasklet)(Matrix* , Matrix* , Matrix* , int, int, int, int, int, 
                        int, int, Sync*, int);

typedef struct tasklet_desc {
  Matrix* t_X;
  Matrix* t_Y;
  Matrix* t_Z;
  int t_x_row;
  int t_x_col;
  int t_y_row;
  int t_y_col;
  int t_z_row;
  int t_z_col;
  int t_n;
} tasklet_desc;

typedef struct work {
    tasklet w_t;
    tasklet_desc w_td;
    Sync* w_sync;
} Work;

class Task { 
public:
    Task(int id, int core) {
        is_done = false;
        is_started = false;
        id_ = id;
        cores_ = core;
        if(core == 1) {
            steal_count_ = core;
        }
        else {
            steal_count_ = ceil(core * log2(core));
        }
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
        /*std::cout << "PUF: " << " x_row: " << w->w_td.t_x_row << " x_col: " 
                    << w->w_td.t_x_col << " y_row: " << w->w_td.t_y_row 
                    << " y_col: " << w->w_td.t_y_col << " z_row: "
                    << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col 
                    << " size: " << w->w_td.t_n<< "\n"; */
        d_.push_front(w);
    }

    void push_back(Work* w) {
        std::unique_lock<std::mutex> lock(m_);
        /*std::cout << "PUB: " << " x_row: " << w->w_td.t_x_row << " x_col: " 
                    << w->w_td.t_x_col << " y_row: " << w->w_td.t_y_row 
                    << " y_col: " << w->w_td.t_y_col << " z_row: " 
                    << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col 
                    << " size: " << w->w_td.t_n << "\n"; */
        d_.push_back(w);
    }

    Work* pop_back() {
        std::unique_lock<std::mutex> lock(m_);
        if(d_.size() == 0) {
            return NULL;
        }
        Work *w = d_.back();
        /*std::cout << "POB: " << " x_row: " << w->w_td.t_x_row << " x_col: " 
                    << w->w_td.t_x_col << " y_row: " << w->w_td.t_y_row 
                    << " y_col: " << w->w_td.t_y_col << " z_row: "
                    << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col 
                    << " size: " << w->w_td.t_n << "\n"; */
       
        d_.pop_back();
        return w;
    }

    Work* pop_front() {
        std::unique_lock<std::mutex> lock(m_);
        if(d_.size() == 0) {
            return NULL;
        }
        Work *w = d_.front();
        /*std::cout << "POF: " << " x_row: " << w->w_td.t_x_row << " x_col: " 
                    << w->w_td.t_x_col << " y_row: " << w->w_td.t_y_row 
                    << " y_col: " << w->w_td.t_y_col << " z_row: " 
                    << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col 
                    << " size: " << w->w_td.t_n << "\n"; */
        d_.pop_front();
        return w;
    }

    int get_deque_size() {
        return d_.size();
    }
    
    void push_sync(Sync* sync) {
       if(sync == NULL) {
           return;
       }
       std::unique_lock<std::mutex> lock(s_m_);
       /*std::cout <<"PUS: "<< "SyncType: " << sync->s_type << " SyncVal: " 
                   << sync->s_val << " x_r: " << sync->s_x_row << " x_c: " 
                   << sync->s_x_col  << " y_r: " << sync->s_y_row << " y_c: " 
                   << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " 
                   << sync->s_z_col << " n: " << sync->s_n << "\n"; */
       st_.push(sync);
    }

    Sync* get_stack_top() {
        std::unique_lock<std::mutex> lock(s_m_);
        if(st_.empty()) {
            return NULL;
        }
        Sync *ret = st_.top();
        /*std::cout <<"ST : "<< "SyncType: " << ret->s_type << " SyncVal: " 
                    << ret->s_val << " x_r: " << ret->s_x_row << " x_c: " 
                    << ret->s_x_col  << " y_r: " << ret->s_y_row << " y_c: " 
                    << ret->s_y_col <<  " z_r: " << ret->s_z_row << " z_c: " 
                    << ret->s_z_col << " n: " << ret->s_n << "\n";*/
        return ret;
    }

    Sync* get_stack_top_ref() {
        std::unique_lock<std::mutex> lock(s_m_);
        if(st_.empty()) {
            return NULL;
        }
        Sync *ret =st_.top();
        /*std::cout <<"STR: "<< "SyncType: " << ret->s_type << " SyncVal: " 
                    << ret->s_val << " x_r: " << ret->s_x_row << " x_c: " 
                    << ret->s_x_col  << " y_r: " << ret->s_y_row << " y_c: " 
                    << ret->s_y_col <<  " z_r: " << ret->s_z_row << " z_c: " 
                    << ret->s_z_col << " n: " << ret->s_n << "\n";*/
        return ret;
    }

    void dec_sync_ref_count(Sync *sync) {
        std::unique_lock<std::mutex> lock(s_m_);
         if(sync == NULL) {
            return;
        }
        /*std::cout <<"DSR: "<< "SyncType: " << sync->s_type << " SyncVal: " 
                    << sync->s_val << " x_r: " << sync->s_x_row << " x_c: " 
                    << sync->s_x_col  << " y_r: " << sync->s_y_row << " y_c: " 
                    << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " 
                    << sync->s_z_col << " n: " << sync->s_n << "\n";*/
        --(sync->s_val); 
        return;
    }

    bool is_sync_ref_count_zero(Sync *sync) {
         std::unique_lock<std::mutex> lock(s_m_);
        if(sync == NULL) {
            return false;
        }
        /*std::cout <<"ISR: "<< "SyncType: " << sync->s_type << " SyncVal: " 
                    << sync->s_val << " x_r: " << sync->s_x_row << " x_c: " 
                    << sync->s_x_col  << " y_r: " << sync->s_y_row << " y_c: " 
                    << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " 
                    << sync->s_z_col << " n: " << sync->s_n << "\n";*/
        return (sync->s_val == 0) ? true : false; 
    }

    int get_sync_ref_count(Sync *sync) {
        std::unique_lock<std::mutex> lock(s_m_);
         if(sync == NULL) {
            return 0;
        }
        /*std::cout <<"GSR: "<< "SyncType: " << sync->s_type << " SyncVal: " 
                    << sync->s_val << " x_r: " << sync->s_x_row << " x_c: " 
                    << sync->s_x_col  << " y_r: " << sync->s_y_row << " y_c: " 
                    << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " 
                    << sync->s_z_col << " n: " << sync->s_n << "\n";*/
        return sync->s_val;
    }

    void pop_sync() {
       std::unique_lock<std::mutex> lock(s_m_);
       if(st_.empty()) {
           return;
       }
       Sync *top = st_.top();
       st_.pop();
    }
    
    bool is_stack_empty() {
        return st_.size() == 0 ? true : false;
    }

    void run() {
        int steal_count = 0;
        while(1) {
            if(!d_.empty()) {
                is_started = true;
                //std::cout << "Inside Run : Popping Job" << "\n";
                Work *w = this->pop_back();
                if(w) {
                    (*(w->w_t))(w->w_td.t_X, w->w_td.t_Y, w->w_td.t_Z, 
                                w->w_td.t_x_row, w->w_td.t_x_col, w->w_td.t_y_row, 
                                w->w_td.t_y_col, w->w_td.t_z_row, w->w_td.t_z_col, 
                                w->w_td.t_n, w->w_sync, id_);
                }
                steal_count = 0;
            }
            else {
                // STEAL
                ++steal_count;
                if(steal_count == 1000*steal_count_) {
                    return;
                }
            }
        }
        return;
    }

private:
    std::mutex m_;
    std::deque<Work *> d_;
    std::stack<Sync *> st_;
    std::mutex s_m_;
    bool is_started;
    bool is_done;
    bool is_top_half_done;
    bool is_bottom_half_done;
    int id_;
    int cores_;
    uint64_t steal_count_; 
};

#endif

