#include <iostream>
#include <math.h>
#include <boost/foreach.hpp>

#include "icetray/I3TrayHeaders.h"
#include "icetray/I3Units.h"
#include "icetray/I3Bool.h"
#include "dataclasses/physics/I3TriggerHierarchy.h"
#include "dataclasses/status/I3TriggerStatus.h"
#include "dataclasses/status/I3DetectorStatus.h"
#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/physics/I3DOMLaunch.h"
#include "dataclasses/TriggerKey.h"
#include "dataclasses/I3Constants.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include "trigger-sim/utilities/DOMSetFunctions.h"
#include "trigger-sim/algorithms/TriggerContainer.h"
#include "trigger-sim/algorithms/TriggerHit.h"
#include "trigger-sim/algorithms/SlowMonopoleTriggerAlgorithm.h"
#include "trigger-sim/utilities/DetectorStatusUtils.h"

/*
 * slow monopole trigger
 */
using namespace std;
using namespace I3Units;

using DetectorStatusUtils::tk_ts_pair_t;
using DetectorStatusUtils::_sourceID;
using DetectorStatusUtils::_typeID;
using DetectorStatusUtils::_configID;
using DetectorStatusUtils::GetTriggerStatus;

const TriggerKey::SourceID SOURCEID(TriggerKey::IN_ICE);
const TriggerKey::TypeID TYPEID(TriggerKey::SLOW_PARTICLE);

SlowMonopoleTriggerAlgorithm::SlowMonopoleTriggerAlgorithm(double t_proximity, double t_min, double t_max, 
                                                           boost::optional<double> deltad, 
                                                           boost::optional<double> alpha_min,
                                                           boost::optional<bool> dc_algo,
                                                           double relv, 
                                                           int min_tuples,
                                                           double max_event_length,
                                                           I3GeometryConstPtr Geometry,
                                                           int domSet, I3MapKeyVectorIntConstPtr customDomSets):
  TriggerService(domSet, customDomSets),
  t_proximity_(t_proximity),     // 2.5 microseconds
  t_min_(t_min),           // 0 microseconds
  t_max_(t_max),           // 500 microseconds
  deltad_(deltad),          // 100 meters  (not used if dc_algo = true)
  alpha_min_(alpha_min),       // 140 deg     (not used if dc_algo = false)
  dc_algo_(dc_algo),           // true since IC2012, false for IC2011
  relv_(relv),            // 0.5 (all values are per definition smaller than 3.0)
  min_tuples_(min_tuples),       // take all by default. pole settings are: 5 for IC2012, 3 for IC2011
  max_event_length_(max_event_length), // 5000 microseconds
  geo_(Geometry)
{
  if(alpha_min_)
    cos_alpha_min_ = cos(alpha_min_.get()*I3Units::degree);

  // just a cross check
  if(dc_algo_ && dc_algo_.get() == true && !deltad_){
    log_error("Misconfigured I3TriggerStatus.");
    log_fatal("If dc_algo is true then delta_d needs to be defined.");
  }

  if(dc_algo_ && dc_algo_.get() == false && !alpha_min_){
    log_error("Misconfigured I3TriggerStatus.");
    log_fatal("If dc_algo is false then alpha_min needs to be defined.");
  }     
}

void SlowMonopoleTriggerAlgorithm::Trigger(){
  // clear global variables before starting
  trigger_container_vector.clear();
  trigger_list.clear();

  // per frame variables
  TriggerHitVector one_hit_list;
  TriggerHitVector two_hit_list;
  double muon_time_window = -1;

  TriggerHitVectorPtr hits(new TriggerHitVector);
  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////   calling trigger    /////////////////////
  ////////////////////////////////////////////////////////////////////////
  for(int i = 0; i < int(hits_->size()); i++)
  {
	TriggerHit payload = hits_->at(i);
	RunTrigger(&one_hit_list, &two_hit_list, &muon_time_window, payload, geo_);
  }

  // Final call of CheckTriggerStatus. This is only necessary in simulation.
  // In the DAQ module the data stream never ends, so the timing condition is fulfilled.
  CheckTriggerStatus(&two_hit_list, geo_);
  log_debug("Trigger status checked. trigger_list.size(): %d", int(trigger_list.size()) );
  return;
}


void SlowMonopoleTriggerAlgorithm::RunTrigger(TriggerHitVector *one_hit_list__, 
				     TriggerHitVector *two_hit_list__, 
				     double *muon_time_window__, 
				     TriggerHit new_hit, 
				     const I3GeometryConstPtr &geo)
{
  if(one_hit_list__->size() == 0) // size is 0, so just add it to the list
    {
      one_hit_list__->push_back(new_hit);
    }
  else // not zero, so compare payload with current one_hit_list if any hlc pair can be formed
    {
      while( fabs(new_hit.time - one_hit_list__->front().time) > 1000.0)
	// makes no sense to compare HLC hits that are 
	//longer apart than 1000 nanoseconds, so remove first from list
        {
            one_hit_list__->erase(one_hit_list__->begin());

            if(one_hit_list__->size() == 0)
            {
                break;
            }
        }

        for(TriggerHitVector::iterator one_hit_iter = one_hit_list__->begin(); 
	    one_hit_iter != one_hit_list__->end(); )
        {
        	// HLC Pair Check
            if(HLCPairCheck(*one_hit_iter, new_hit))  
	      // the geometry here was just a test and can be deleted later
            {
	      // the payload is the first hit from the HLC pair
            	TriggerHit check_payload = *one_hit_iter;  
                if(two_hit_list__->size() == 0)        // the pair list is empty
                {
                    if(*muon_time_window__ == -1)
                    {
                        two_hit_list__->push_back(check_payload);
                    }
                    else
                    {
                        if(check_payload.time - *muon_time_window__ <= t_proximity_)
                        {
                            *muon_time_window__ = check_payload.time;
                        }
                        else
                        {
                            two_hit_list__->push_back(check_payload);
                            *muon_time_window__ = -1;
                        }
                    }
                }
                else // the pair list is not empty
                {
                    if(*muon_time_window__ == -1)
                    {
                        if(check_payload.time - two_hit_list__->back().time <= t_proximity_)
                        {
                            *muon_time_window__ = check_payload.time;
                            two_hit_list__->pop_back();
                        }
                        else
                        {
                            if((check_payload.time - two_hit_list__->back().time < t_max_)
                               && (check_payload.time - two_hit_list__->front().time 
				   < max_event_length_))
                            {
                                two_hit_list__->push_back(check_payload);
                            }
                            else
                            {
			      // checks current two_hit_list for 3-tuples
                                CheckTriggerStatus(two_hit_list__, geo); 
                                two_hit_list__->push_back(check_payload);
                            }
                        }
                    }
                    else
                    {
                        if(check_payload.time - *muon_time_window__ <= t_proximity_)
                        {
                            *muon_time_window__ = check_payload.time;
                        }
                        else
                        {
                            *muon_time_window__ = -1;
                            if((check_payload.time - two_hit_list__->back().time < t_max_)
                               && (check_payload.time - two_hit_list__->front().time < max_event_length_))
                            {
                                two_hit_list__->push_back(check_payload);
                            }
                            else
                            {
                                CheckTriggerStatus(two_hit_list__, geo);
                                two_hit_list__->push_back(check_payload);
                            }
                        }
                    }
                }

		// deletes the hit and increments the iterator
                one_hit_iter = one_hit_list__->erase(one_hit_iter); 	
            }
            else  // no HLC found
            {
	      // just increment the iterator
                one_hit_iter++ ;
            }
        }

	// at the end add the current hitPayload for further comparisons
        one_hit_list__->push_back(new_hit); 

        if(two_hit_list__->size() > 0) // definitely cannot produce a trigger
        {
        	if( (one_hit_list__->front().time - two_hit_list__->back().time) > t_max_)
        	{
                CheckTriggerStatus(two_hit_list__, geo);
        	}
        }
    }
}

bool SlowMonopoleTriggerAlgorithm::HLCPairCheck(TriggerHit hit1, TriggerHit hit2)
{
	int string_nr1 = hit1.string;
	int string_nr2 = hit2.string;

    if(string_nr1 == string_nr2)
    {
        int om_nr1 = hit1.pos;
        int om_nr2 = hit2.pos;

        if( abs(om_nr1 - om_nr2) <= 2)
        {
            return true;
        }
    }
    return false;
}

void SlowMonopoleTriggerAlgorithm::CheckTriggerStatus(TriggerHitVector *two_hit_list__, 
					     const I3GeometryConstPtr &geo)
{
    int list_size = two_hit_list__->size();
    if(list_size >= 3)
    {
        for(TriggerHitVector::iterator iter_1 = two_hit_list__->begin(); 
	    iter_1 != two_hit_list__->end() -2 ; iter_1++ )
        {
        	for(TriggerHitVector::iterator iter_2 = iter_1 + 1; 
		    iter_2 != two_hit_list__->end() - 1 ; iter_2++ )
        	{
        		for(TriggerHitVector::iterator iter_3 = iter_2 + 1 ; 
			    iter_3 != two_hit_list__->end(); iter_3++ )
        		{
        			CheckTriple(*iter_1, *iter_2, *iter_3, geo);
        		}
        	}
        }
    }


    for(TriggerContainerVector::iterator trg_iter = trigger_container_vector.begin(); 
	trg_iter != trigger_container_vector.end(); trg_iter ++)
    {
    	TriggerContainer slowmptrigger = *trg_iter;

        // Form trigger with each slowmptrigger object, if num_tuples is fulfilled
        if(slowmptrigger.GetNTuples() >= min_tuples_)
        {
            log_debug("FOUND TRIGGER: start: %f, length :%f with %d tuples.",
		     slowmptrigger.GetTrigger()->GetTriggerTime(), 
		     slowmptrigger.GetTrigger()->GetTriggerLength(), 
		     slowmptrigger.GetNTuples() );
            trigger_list.push_back(*slowmptrigger.GetTrigger());
            
        }
    }
    // Convert these all into the same format as the other trigger algorithms.
    BOOST_FOREACH(TriggerContainer container, trigger_container_vector){
        I3TriggerPtr trigger = container.GetTrigger();

        log_warn("SLOP trigger!");

        // Start time defined by a hit.
        TriggerHit start(0,0,0,1);
        TriggerHit end(0,0,0,1);

        start.time = trigger->GetTriggerTime();
        end.time = start.time + trigger->GetTriggerLength();

        // Add them to the list of triggers
        TriggerHitVector trig_vec;
        trig_vec.push_back(start);
        trig_vec.push_back(end);
        triggers_.push_back(trig_vec);
        triggerCount_++;
    }
    trigger_container_vector.clear();
    two_hit_list__->clear();
}

void SlowMonopoleTriggerAlgorithm::CheckTriple(TriggerHit hit1, 
				      TriggerHit hit2,  
				      TriggerHit hit3, 
				      const I3GeometryConstPtr &geo)
{
  double t_diff1 = hit2.time - hit1.time;
  double t_diff2 = hit3.time - hit2.time;
  if((t_diff1 > t_min_) && (t_diff2 > t_min_) && (t_diff1 < t_max_) && (t_diff2 < t_max_))
    {
      double t_diff3 = hit3.time - hit1.time;
      
      double p_diff1 = getDistance(hit1, hit2, geo);
      double p_diff2 = getDistance(hit2, hit3, geo);
      double p_diff3 = getDistance(hit1, hit3, geo);
      log_debug("    ->step2 - p_diff1: %f, p_diff2: %f, p_diff3: %f", 
		p_diff1, p_diff2, p_diff3);
      
      if ( !( (p_diff1 > 0) && (p_diff2 > 0) && (p_diff3 > 0) ))
        {
	  log_debug("exiting check triple because p_diff1: %f, p_diff2: %f, p_diff3: %f", 
		    p_diff1, p_diff2, p_diff3);
	  return;
        }
      
      double cos_alpha = ( pow(p_diff1,2) + pow(p_diff2,2) - pow(p_diff3,2) ) 
	/ ( 2*p_diff1*p_diff2 );

      // being lazy...
      if(alpha_min_){
	log_debug("alpha_min_: %f deg (=%f), cos_alpha_min_: %f,"
		  " alpha: %f deg (=%f), cos_alpha: %f", 
		  alpha_min_.get(), alpha_min_.get()*I3Units::degree, cos_alpha_min_, 
		  acos(cos_alpha)/I3Units::degree, acos(cos_alpha), cos_alpha);
      }
      if(   ( deltad_ && (p_diff1 + p_diff2 - p_diff3 <= deltad_.get()) )   
	    ||   ( alpha_min_ && (cos_alpha <= cos_alpha_min_) )   )
        {
	  double inv_v1 = t_diff1 / p_diff1;
	  double inv_v2 = t_diff2 / p_diff2;
	  double inv_v3 = t_diff3 / p_diff3;
	  
	  log_debug("dc_algo_ is set to %d cos_alpha is %f ", dc_algo_.get() , cos_alpha);
	  
	  double inv_v_mean = (inv_v1 + inv_v2 + inv_v3)/3.0;
	  
	  log_debug("inv_v1: %f, inv_v2: %f, inv_v3: %f, inv_v_mean: %f", 
		    inv_v1, inv_v2, inv_v3, inv_v_mean);
	  
	  log_debug("        ->step3 - inv_v_mean %f: " , fabs(inv_v2 - inv_v1) / inv_v_mean);
	  if(fabs(inv_v2 - inv_v1) / inv_v_mean <= relv_)
            {
	      // Found Triple
	      log_debug("Found Triple: t1: %f, t2: %f, t3: %f, t_diff1: %f  t_diff2: %f", 
			hit1.time, hit2.time, hit3.time, t_diff1, t_diff2);
	      
	      double triple_start = hit1.time;
	      double triple_end = hit3.time;

	      //  prepare a lot of stuff which is going to be written with the trigger	      
	      if(trigger_container_vector.size() == 0)
                {
		  I3TriggerPtr new_trig = I3TriggerPtr(new I3Trigger);
		  new_trig->SetTriggerFired(true);
		  new_trig->SetTriggerTime(triple_start);
		  new_trig->SetTriggerLength(triple_end-triple_start);
		  
		  TriggerContainer new_cont(new_trig);
		  trigger_container_vector.push_back(new_cont);
		  
		  log_debug("** Adding first trigger to list **");
                }
	      else
                {
		  double trigger_start_temp = 
		    trigger_container_vector.back().GetTrigger()->GetTriggerTime();

		  double trigger_end_temp = 
		    trigger_start_temp 
		    + trigger_container_vector.back().GetTrigger()->GetTriggerLength();
		  
		  log_debug("Trigger TEMP: %f %f", trigger_start_temp, trigger_end_temp);
		  
		  if((triple_start >= trigger_start_temp) 
		     && (triple_start <= trigger_end_temp) 
		     && (triple_end > trigger_end_temp)) //overlap
                    {
		      trigger_container_vector.back().GetTrigger()->
			SetTriggerLength(triple_end-trigger_start_temp);

		      trigger_container_vector.back().IncreaseNTuples();
		      log_debug("** Found overlap with existing trigger **");
                    }
		  else if((triple_start >= trigger_start_temp) 
			  && (triple_end <= trigger_end_temp)) // contained tuple
                    {
		      trigger_container_vector.back().IncreaseNTuples();
		      log_debug("** Found contained tuple **");

                    }
		  else if(triple_start > trigger_end_temp)
                    {
		      
		      I3TriggerPtr new_trig = I3TriggerPtr(new I3Trigger);
		      new_trig->SetTriggerFired(true);
		      new_trig->SetTriggerTime(triple_start);
		      new_trig->SetTriggerLength(triple_end-triple_start);
		      
		      TriggerContainer new_cont(new_trig);
		      
		      trigger_container_vector.push_back(new_cont);
		      
		      log_debug("** Found SECOND TRIGGER **");
                    }
                }
            }
        }
    }
}

double SlowMonopoleTriggerAlgorithm::getDistance(TriggerHit hit1, 
					TriggerHit hit2, 
					const I3GeometryConstPtr &geo)
{
  I3OMGeoMap::const_iterator geo_iterator_1 = geo->omgeo.find(OMKey(hit1.string, hit1.pos));
  I3OMGeoMap::const_iterator geo_iterator_2 = geo->omgeo.find(OMKey(hit2.string, hit2.pos));
  
  double x1 = geo_iterator_1->second.position.GetX();
  double y1 = geo_iterator_1->second.position.GetY();
  double z1 = geo_iterator_1->second.position.GetZ();
  double x2 = geo_iterator_2->second.position.GetX();
  double y2 = geo_iterator_2->second.position.GetY();
  double z2 = geo_iterator_2->second.position.GetZ();
  
  double diff = sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2) );
  
  return diff;
}