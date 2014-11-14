#ifndef __ncm_qapt_server_h__
#define __ncm_qapt_server_h__

//:
// \file
// \brief Qt interface to an APT motor controller server class.
//
// \author Phil Tresadern

#include <vcl_iostream.h>
#include <vcl_map.h>

#include <QObject>
#include <QTimer>

#include <nailfold/ncm_apt_server.h>

class ncm_qapt_server : public QObject, 
                        public ncm_apt_server
{
  Q_OBJECT // needed if we want to use signals or slots

//  INTERFACE

public:

  enum apt_axis {
    AptInvalidFirst = -1,
    AptX = 0,
    AptY,
    AptZ,

    First = AptX,
    Last = AptZ
  }; 

  //: Default constructor.
  ncm_qapt_server();


signals:

  void move_complete(ncm_qapt_server::apt_axis);

public slots:

  //: Set the datum point for limited movements.
  void set_reference_position(ncm_qapt_server::apt_axis axis);

  //: Set the maximum displacement for limited movements.
  void set_max_displacement(ncm_qapt_server::apt_axis axis, 
                            float max_displacement = 1e9);

	//: Populate controllers_ with new instances of ncm_apt_controller.
  void create_controllers();

	//: Clear controllers_.
  void destroy_controllers();


//  IMPLEMENTATION

private slots:

  // The following slots are defined for commands that do not return 
  // immediately and that may hang the system (i.e. those calling 
  // MoveRelative*, MoveAbsolute* and MoveHome).

  // Because they should be called using signals, rather than directly,
  // I've made them private.

  //: Move specified axis to the zero position (with respect to the reverse 
  //  limit switch).
  void move_zero(ncm_qapt_server::apt_axis axis, 
                 bool return_immediately = true);

  //: Move specified axis to the home position.
  void move_home(ncm_qapt_server::apt_axis axis, 
                 float home_position, 
                 bool return_immediately = true);

  //: Move specified axis to the home position.
  void move_home(ncm_qapt_server::apt_axis axis, 
                 bool return_immediately = true);

  //: Move specified axis to given position (relative to home position).
  void move_to(ncm_qapt_server::apt_axis axis, 
               float position,
               bool return_immediately = true);

  //: Move specified axis by a given distance.
  void move_by(ncm_qapt_server::apt_axis axis, 
               float distance,
               bool return_immediately = true);

  //: Monitor displacement if need be.
  void on_apt_timer();


private:
  
  //: Return a pointer to an axis, given the name.
  ncm_apt_controller_base* get_axis(ncm_qapt_server::apt_axis axis);


  //: Timer that monitors position and emits signals where necessary.
  QTimer apt_timer_;

  //: Reference position from which max_displacement is measured.
  vcl_map<apt_axis, float> reference_position_;

  //: Upper limit of motion during a move_at_velocity() command.
  vcl_map<apt_axis, float> max_displacement_;
};

#endif // __ncm_qapt_server_h__
