#include <windows.h>
#include <APTAPI.h>

#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vcl_cmath.h>
#include <vcl_ctime.h>

#include <nailfold/ncm_apt_server.h>
#include <nailfold/ncm_apt_controller.h>

void test_apt_server()
{
  ncm_apt_server apt_server;

  apt_server.set_X_id(83837992);

  // Move to zero
  vcl_cout << "Zeroing...";
  apt_server.X()->move_zero(false);
  vcl_cout << "done" << vcl_endl;

  vcl_cout << "Homing...";
  apt_server.X()->move_home(4.0f, false);
  vcl_cout << "done" << vcl_endl;

  // Move back 1mm
  apt_server.X()->move_to(-1.0f, false);

  float velocity = 0.0f;
  float acceleration = 1.0f;

  vcl_ofstream ofs("u:/tmp/positions.txt");

  int t0 = vcl_clock();

  // Execute simple harmonic motion
  unsigned nSteps = 800;
  while(nSteps != 0)
  {
    // Stop gracefully if actuator gets too close to its limits.
    if (apt_server.X()->distance_to_endstop() < 1.0f)
    {
      apt_server.X()->stop();
      break;
    }

    // Calculate the velocity at which we want to move.
    float dt = 0.05f; // Probably not correct
    float acceleration = -apt_server.X()->position();
    velocity = velocity + acceleration*dt; // v = u + at

    // Limit the velocity (since our estimates are a bit crap).
    if (velocity > 1.0f)
      velocity = 1.0;
    else if (velocity < -1.0f)
      velocity = -1.0;

    apt_server.X()->move_at_velocity(velocity);

    int t = vcl_clock();
    if (t-t0 > 10)
    {
      ofs << t << "\t" 
          << apt_server.X()->absolute_position() << "\t" 
          << velocity << "\t"
          << acceleration
          << vcl_endl;
      --nSteps;
      t0 = t;
    }
  }

  apt_server.X()->stop();
  ofs.close();
}

void test_apt_limits()
{
  const int id_ = 83837992;
  const float home_velocity = 1.0f;
  const float home_position = 3.0f;
  const bool bWait = true;

  APTInit();
  EnableEventDlg(false);

  InitHWDevice(id_);
  MOT_SetHomeParams(id_, HOME_REV, HOMELIMSW_REV,
                    home_velocity, home_position);
  MOT_MoveHome(id_, bWait);

  float min_position_, max_position_, position;
  long units;
  float pitch;
  MOT_GetStageAxisInfo(id_, &min_position_, &max_position_, 
                            &units, &pitch);
  min_position_ -= home_position;
  max_position_ -= home_position;
  MOT_SetStageAxisInfo(id_, min_position_, max_position_, 
                            units, pitch);

  MOT_GetPosition(id_, &position);
  vcl_cout << "min: " << min_position_ << "; "
           << "max: " << max_position_ << "; " 
           << "curr: " << position << vcl_endl;

  MOT_MoveRelativeEx(id_, 1.0f, bWait);

  MOT_GetPosition(id_, &position);
  vcl_cout << "min: " << min_position_ << "; "
           << "max: " << max_position_ << "; " 
           << "curr: " << position << vcl_endl;

  const float zero_offset = 0.0f;
  MOT_SetHomeParams(id_, HOME_REV, HOMELIMSW_REV,
                    home_velocity, zero_offset);

  MOT_GetPosition(id_, &position);
  vcl_cout << "min: " << min_position_ << "; "
           << "max: " << max_position_ << "; " 
           << "curr: " << position << vcl_endl;

  APTCleanUp();
}

int main()
{
  //test_apt_limits();

  test_apt_server();

  return 0;
}

