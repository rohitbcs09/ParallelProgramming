#ifndef _MAT_MUL_STEAL_H
#define _MAT_MUL_STEAL_H

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <deque>
#include <stack>
#include <stdint.h>

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

Sync* create_sync_state(Sync *p_sync, int type, int x1, int x2, int y1, int y2, int z1, int z2, int size) {
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
        /*std::cout << "PUF: " << " x_row: " << w->w_td.t_x_row << " x_col: " << w->w_td.t_x_col
                  << " y_row: " << w->w_td.t_y_row << " y_col: " << w->w_td.t_y_col << " z_row: "
                  << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col << " size: " << w->w_td.t_n<< "\n"; */
        d_.push_front(w);
    }

    void push_back(Work* w) {
        std::unique_lock<std::mutex> lock(m_);
        /*std::cout << "PUB: " << " x_row: " << w->w_td.t_x_row << " x_col: " << w->w_td.t_x_col
                  << " y_row: " << w->w_td.t_y_row << " y_col: " << w->w_td.t_y_col << " z_row: "
                  << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col << " size: " << w->w_td.t_n << "\n"; */
       
        d_.push_back(w);
    }

    Work* pop_back() {
        std::unique_lock<std::mutex> lock(m_);
		if(d_.size() == 0) {
            return NULL;
		}
        Work *w = d_.back();
        /*std::cout << "POB: " << " x_row: " << w->w_td.t_x_row << " x_col: " << w->w_td.t_x_col
                  << " y_row: " << w->w_td.t_y_row << " y_col: " << w->w_td.t_y_col << " z_row: "
                  << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col << " size: " << w->w_td.t_n << "\n"; */
       
        d_.pop_back();
        return w;
    }

    Work* pop_front() {
        std::unique_lock<std::mutex> lock(m_);
		if(d_.size() == 0) {
            return NULL;
		}
        Work *w = d_.front();
        /*std::cout << "POF: " << " x_row: " << w->w_td.t_x_row << " x_col: " << w->w_td.t_x_col
                  << " y_row: " << w->w_td.t_y_row << " y_col: " << w->w_td.t_y_col << " z_row: "
                  << w->w_td.t_z_row << " z_col: " << w->w_td.t_z_col << " size: " << w->w_td.t_n << "\n"; */
        d_.pop_front();
        return w;
    }

    int get_deque_size() {
        return d_.size();
    }
    
    void push_sync(Sync* sync) {
       std::unique_lock<std::mutex> lock(s_m_);
	   /*std::cout <<"PUS: "<< "SyncType: " << sync->s_type << " SyncVal: " << sync->s_val 
	             << " x_r: " << sync->s_x_row << " x_c: " << sync->s_x_col  << " y_r: " << sync->s_y_row 
				 << " y_c: " << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " << sync->s_z_col 
				 << " n: " << sync->s_n << "\n"; */
       st_.push(sync);
    }

    Sync stack_top() {
        std::unique_lock<std::mutex> lock(s_m_);
        Sync *ret = st_.top();
	    /*std::cout <<"ST : "<< "SyncType: " << ret->s_type << " SyncVal: " << ret->s_val
	             << " x_r: " << ret->s_x_row << " x_c: " << ret->s_x_col  << " y_r: " << ret->s_y_row 
				 << " y_c: " << ret->s_y_col <<  " z_r: " << ret->s_z_row << " z_c: " << ret->s_z_col 
				 << " n: " << ret->s_n << "\n";*/
        return *ret;
    }

    Sync* stack_top_ref() {
        std::unique_lock<std::mutex> lock(s_m_);
        Sync *ret =st_.top();
	    /*std::cout <<"STR: "<< "SyncType: " << ret->s_type << " SyncVal: " << ret->s_val
                 << " x_r: " << ret->s_x_row << " x_c: " << ret->s_x_col  << " y_r: " << ret->s_y_row 
			   	 << " y_c: " << ret->s_y_col <<  " z_r: " << ret->s_z_row << " z_c: " << ret->s_z_col 
				 << " n: " << ret->s_n << "\n";*/

        return ret;
    }

    Sync* dec_sync_ref(Sync *sync) {
        std::unique_lock<std::mutex> lock(s_m_);
 	    if(sync == NULL) {
		    return NULL;
		}
	    /*std::cout <<"DSR: "<< "SyncType: " << sync->s_type << " SyncVal: " << sync->s_val
	             << " x_r: " << sync->s_x_row << " x_c: " << sync->s_x_col  << " y_r: " << sync->s_y_row 
				 << " y_c: " << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " << sync->s_z_col 
				 << " n: " << sync->s_n << "\n";*/
        (sync->s_val)--; 
		if(sync->s_val == 0) {
	       /*std::cout <<"POS: "<< "SyncType: " << sync->s_type << " SyncVal: " << sync->s_val
	             << " x_r: " << sync->s_x_row << " x_c: " << sync->s_x_col  << " y_r: " << sync->s_y_row 
				 << " y_c: " << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " << sync->s_z_col 
				 << " n: " << sync->s_n << "\n";*/
    
		   st_.pop();
		}
		return sync;
    }

    bool is_sync_ref_zero(Sync *sync) {
 	    std::unique_lock<std::mutex> lock(s_m_);
        if(sync == NULL) {
		    return true;
		}
	    /*std::cout <<"ISR: "<< "SyncType: " << sync->s_type << " SyncVal: " << sync->s_val
	             << " x_r: " << sync->s_x_row << " x_c: " << sync->s_x_col  << " y_r: " << sync->s_y_row 
				 << " y_c: " << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " << sync->s_z_col 
				 << " n: " << sync->s_n << "\n";*/
        return (sync->s_val == 0) ? true : false; 
    }
   
    int get_sync_ref_count(Sync *sync) {
        std::unique_lock<std::mutex> lock(s_m_);
 	    if(sync == NULL) {
		    return 0;
		}
	    /*std::cout <<"GSR: "<< "SyncType: " << sync->s_type << " SyncVal: " << sync->s_val 
	             << " x_r: " << sync->s_x_row << " x_c: " << sync->s_x_col  << " y_r: " << sync->s_y_row 
				 << " y_c: " << sync->s_y_col <<  " z_r: " << sync->s_z_row << " z_c: " << sync->s_z_col 
				 << " n: " << sync->s_n << "\n";*/
 
        return sync->s_val;
    }

    void pop_sync() {
       std::unique_lock<std::mutex> lock(s_m_);
	   Sync *top = st_.top();
	   if(top) {}
	       /*std::cout <<"POS: "<< "SyncType: " << top->s_type << " SyncVal: " << top->s_val 
	             << " x_r: " << top->s_x_row << " x_c: " << top->s_x_col  << " y_r: " << top->s_y_row 
				 << " y_c: " << top->s_y_col <<  " z_r: " << top->s_z_row << " z_c: " << top->s_z_col 
				 << " n: " << top->s_n << "\n";*/
 

	   else {}
	      /* std::cout <<"POS: "<< "SyncType: " << -1 << " SyncVal: " << -1
	             << " x_r: " << top->s_x_row << " x_c: " << top->s_x_col  << " y_r: " << top->s_y_row 
				 << " y_c: " << top->s_y_col <<  " z_r: " << top->s_z_row << " z_c: " << top->s_z_col 
				 << " n: " << top->s_n << "\n"; */
 

       st_.pop();
    }
    
    bool is_stack_empty() {
        return st_.size() == 0 ? true : false;
    }

    Work* steal_random_work();

    void run() {
        int count = 0;
        while(1) {
            if(!d_.empty()) {
			    is_started = true;
                //std::cout << "Inside Run : Popping Job" << "\n";
                Work *w = this->pop_back();
                (*(w->w_t))(w->w_td.t_X, w->w_td.t_Y, w->w_td.t_Z, w->w_td.t_x_row,
                          w->w_td.t_x_col, w->w_td.t_y_row, w->w_td.t_y_col, w->w_td.t_z_row, 
                          w->w_td.t_z_col, w->w_td.t_n, w->w_sync, id_);
            }
            else {
                Work *w = steal_random_work();
				if(w != NULL) {
                    (*(w->w_t))(w->w_td.t_X, w->w_td.t_Y, w->w_td.t_Z, w->w_td.t_x_row,
                               w->w_td.t_x_col, w->w_td.t_y_row, w->w_td.t_y_col, w->w_td.t_z_row, 
                               w->w_td.t_z_col, w->w_td.t_n, w->w_sync, id_);
	                count = 0;
			    }
				else {
                    ++count;
                    if(count == 100) {
                        std::cout << "Inside Run : Steal Failure after 20 attempts" << "\n";
                        return;
				    }
                }
                // STEAL
            }
        }
        return;
    }

private:
    std::deque<Work *> d_;
    std::stack<Sync *> st_;
    std::mutex m_;
    std::mutex s_m_;
    bool is_started;
    bool is_done;
    int id_;
};

#endif

