#include "time_series_dd.h"

namespace shyft::time_series::dd {
    
        shared_ptr<ipoint_ts> abin_op_scalar_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;// do not propagate ref_counting more than once
                return c.evaluate(rhs.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_rhs= c.evaluate(rhs.ts);
            abin_op_scalar_ts tmp;
            tmp.lhs=lhs;tmp.op=op;tmp.rhs=apoint_ts(eval_rhs);tmp.ta=ta;tmp.fx_policy=fx_policy;tmp.bound=bound;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),fx_policy);
            c.register_ts(this,r);
            return r;
        }
        
        shared_ptr<ipoint_ts> abin_op_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
             if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 c.evaluate(lhs.ts);
                return c.evaluate(rhs.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_rhs= c.evaluate(rhs.ts);
            auto eval_lhs =c.evaluate(lhs.ts);
            abin_op_ts tmp;
            tmp.lhs=apoint_ts(eval_lhs);tmp.op=op;tmp.rhs=apoint_ts(eval_rhs);tmp.ta=ta;tmp.fx_policy=fx_policy;tmp.bound=bound;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),fx_policy);
            c.register_ts(this,r);
            return r;
        }

        shared_ptr<ipoint_ts> abin_op_ts_scalar::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
             if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                return c.evaluate(lhs.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_lhs =c.evaluate(lhs.ts);
            abin_op_ts_scalar tmp;
            tmp.lhs=apoint_ts(eval_lhs);tmp.op=op;tmp.rhs=rhs;tmp.ta=ta;tmp.fx_policy=fx_policy;tmp.bound=bound;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),fx_policy);
            c.register_ts(this,r);
            return r;
        }

        shared_ptr<ipoint_ts> abs_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
             if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            abs_ts tmp;
            tmp.ts=eval_ts;tmp.ta=ta;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),ts->point_interpretation());
            c.register_ts(this,r);
            return r;
        }
        
        shared_ptr<ipoint_ts> accumulate_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
             if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                return c.evaluate(ts);
            }
            
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            accumulate_ts tmp;
            tmp.ts=eval_ts;tmp.ta=ta;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }
         
       shared_ptr<ipoint_ts> aref_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            return rep;
        }

       shared_ptr<ipoint_ts> average_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            average_ts tmp;
            tmp.ts=eval_ts;tmp.ta=ta;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }
        
        shared_ptr<ipoint_ts> bucket_ts::evaluate(eval_ctx& c, shared_ptr<ipoint_ts> const&) const {
            if (auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                return c.evaluate(ts);
            }
            if (c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts = c.evaluate(ts);
            bucket_ts tmp;
            tmp.ts = eval_ts; tmp.ta = ta;tmp.p=p; tmp.bound = bound;
            auto r = make_shared<gpoint_ts>(ta, tmp.values(), point_interpretation());
            c.register_ts(this, r);
            return r;
        }


       shared_ptr<ipoint_ts> decode_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            decode_ts tmp;
            tmp.ts=eval_ts;tmp.p=p;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }

        shared_ptr<ipoint_ts> convolve_w_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts_impl.ts.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts_impl.ts.ts);
            convolve_w_ts tmp;
            tmp.ts_impl.ts=apoint_ts(eval_ts);tmp.ts_impl.w=ts_impl.w;tmp.ts_impl.policy=ts_impl.policy;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }

       shared_ptr<ipoint_ts> derivative_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            derivative_ts tmp;
            tmp.ts=eval_ts;tmp.dm=dm;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }
        
        shared_ptr<ipoint_ts> extend_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                c.evaluate(rhs.ts);
                 return c.evaluate(lhs.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_rhs= c.evaluate(rhs.ts);
            auto eval_lhs =c.evaluate(lhs.ts);
            extend_ts tmp;
            tmp.lhs=apoint_ts(eval_lhs);
            tmp.rhs=apoint_ts(eval_rhs);
            tmp.ets_split_p=ets_split_p;
            tmp.split_at=split_at;
            tmp.ets_fill_p=ets_fill_p;
            tmp.fill_value=fill_value;
            tmp.ta=ta;
            tmp.fx_policy=fx_policy;
            tmp.bound=bound;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),fx_policy);
            c.register_ts(this,r);
            return r;
        }
        
        shared_ptr<ipoint_ts> gpoint_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            return shared_this;// trival case, evaluate to shared_this always
        }
        
        shared_ptr<ipoint_ts> ice_packing_recession_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                c.evaluate(flow_ts.ts);
                 return c.evaluate(ice_packing_ts.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_flow_ts=c.evaluate(flow_ts.ts);
            auto eval_ice_packing_ts =c.evaluate(ice_packing_ts.ts);
            ice_packing_recession_ts tmp;
            tmp.flow_ts=apoint_ts(eval_flow_ts);
            tmp.ice_packing_ts=apoint_ts(eval_ice_packing_ts);
            tmp.fx_policy=fx_policy;
            tmp.bound=bound;
            tmp.ipr_param=ipr_param;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),fx_policy);
            c.register_ts(this,r);
            return r;
        }

        shared_ptr<ipoint_ts> ice_packing_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts.temp_ts.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts=c.evaluate(ts.temp_ts.ts);
            ice_packing_ts tmp;
            tmp.ts.temp_ts=apoint_ts(eval_ts);
            tmp.ts.ip_param=ts.ip_param;
            tmp.ts.ipt_policy=ts.ipt_policy;
            tmp.ts.bound=ts.bound;
            tmp.ts.fx_policy=ts.fx_policy;
            
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }

       shared_ptr<ipoint_ts> inside_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            inside_ts tmp;
            tmp.ts=eval_ts;tmp.p=p;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }

       shared_ptr<ipoint_ts> integral_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            integral_ts tmp;
            tmp.ts=eval_ts;tmp.ta=ta;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }
        
       shared_ptr<ipoint_ts> krls_interpolation_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts.ts);
            krls_interpolation_ts tmp;
            tmp.ts=apoint_ts(eval_ts);tmp.predictor=predictor;tmp.bound=bound;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }

       shared_ptr<ipoint_ts> periodic_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return nullptr;
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto r= make_shared<gpoint_ts>(time_axis(),ts.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }

       shared_ptr<ipoint_ts> qac_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 c.evaluate(cts);
                 return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            shared_ptr<ipoint_ts> eval_cts= cts?c.evaluate(cts):cts;
            qac_ts tmp;
            tmp.ts=eval_ts;
            tmp.cts=cts;
            tmp.p=p;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }

       shared_ptr<ipoint_ts> rating_curve_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts.level_ts.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts.level_ts.ts);
            rating_curve_ts tmp;
            tmp.ts.level_ts=apoint_ts(eval_ts);
            tmp.ts.rc_param=ts.rc_param;
            tmp.ts.fx_policy=ts.fx_policy;
            tmp.ts.bound=ts.bound;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }
        
       shared_ptr<ipoint_ts> time_shift_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                 return c.evaluate(ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_ts =c.evaluate(ts);
            time_shift_ts tmp;
            tmp.ts=eval_ts;
            tmp.dt=dt;
            tmp.ta=ta;
            auto r= make_shared<gpoint_ts>(time_axis(),tmp.values(),point_interpretation());
            c.register_ts(this,r);
            return r;
        }
        
        shared_ptr<ipoint_ts> use_time_axis_from_ts::evaluate(eval_ctx&c, shared_ptr<ipoint_ts> const& shared_this) const {
            if(auto ref_c=c.ref_counting(this)) { if(ref_c>1u) return nullptr;
                c.evaluate(rhs.ts);
                 return c.evaluate(lhs.ts);
            }
            if(c.is_evaluated(this))
                return c.evaluated[this];
            auto eval_rhs= c.evaluate(rhs.ts);
            auto eval_lhs =c.evaluate(lhs.ts);
            use_time_axis_from_ts tmp;
            tmp.lhs=apoint_ts(eval_lhs);
            tmp.rhs=apoint_ts(eval_rhs);//strictly not needed, we are only for the time-axis, that should be in place by now
            tmp.ta=ta;
            tmp.fx_policy=fx_policy;
            tmp.bound=bound;
            auto r= make_shared<gpoint_ts>(ta,tmp.values(),fx_policy);
            c.register_ts(this,r);
            return r;
        }        
}
