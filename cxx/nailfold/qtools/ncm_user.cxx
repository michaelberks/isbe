#include "ncm_user.h"

//: Define class member functions

//: Constructors
ncm_user::ncm_user(int id /*=1*/)
: name_("default"),
	host_(""),
  user_id_(id) //
{
}

ncm_user::~ncm_user()
{
}

//
//: Return user name
QString ncm_user::name() const
{
  return name_;
}

//
//: Return user id
int ncm_user::user_id() const
{
  return user_id_;
}

//
//: Set user name
void ncm_user::set_name(QString name)
{
  name_ = name;
}

//
//: Set host
void ncm_user::set_host(QString host)
{
  host_ = host;
}
