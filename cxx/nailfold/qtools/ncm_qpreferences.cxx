#include "ncm_qpreferences.h"

#include <nailfold/ncm_vessel.h>

//
//  Define class member functions

//: Constructor
ncm_qpreferences::ncm_qpreferences()
: rezoom_factor_(-1)
{
}

//  Apex labelling parameters

//
//: By how much to rezoom when centring on a vessel (-1 = don't rezoom)
void ncm_qpreferences::set_rezoom_factor(int rezoom_factor)
{
  rezoom_factor_ = rezoom_factor;
}
int ncm_qpreferences::rezoom_factor() const
{
  return rezoom_factor_;
}

//
//: Maximum and minimum distance between consecutive points on a vessel path
void ncm_qpreferences::set_vessel_minimum_interpoint_distance(double distance)
{
  ncm_vessel::set_minimum_inter_point_distance(distance);
}
double ncm_qpreferences::vessel_minimum_interpoint_distance() const
{
  return ncm_vessel::minimum_inter_point_distance();
}
void ncm_qpreferences::set_vessel_maximum_interpoint_distance(double distance)
{
  ncm_vessel::set_maximum_inter_point_distance(distance);
}
double ncm_qpreferences::vessel_maximum_interpoint_distance() const
{
  return ncm_vessel::maximum_inter_point_distance();
}