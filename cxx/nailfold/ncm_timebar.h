#ifndef __NCM_TIMEBAR_H__
#define __NCM_TIMEBAR_H__

class ncm_timebar
{
public: // methods

  ncm_timebar(vcl_string title, unsigned n_ticks);

  void reset(vcl_string title, unsigned n_ticks);

  void advance(unsigned n = 1);

private: // variables

  unsigned n_ticks_;
  unsigned n_chars_;

  unsigned progress_;
};

inline
ncm_timebar::ncm_timebar(vcl_string title, unsigned n_ticks)
{
  reset(title, n_ticks);
}

inline
void ncm_timebar::reset(vcl_string title, unsigned n_ticks)
{
  title += " ";
  vcl_cout << title;

  n_ticks_ = n_ticks;

  n_chars_ =   80 // number of characters in a line
             - title.length() // #chars used by the title
             - 1; // leave a space at the end of the line

  progress_ = 0;
}

inline
void ncm_timebar::advance(unsigned n /* = 1 */)
{
  double progress0 = n_chars_ * static_cast<double>(progress_) / n_ticks_;
  
  progress_ += n;
  double progress1 = n_chars_ * static_cast<double>(progress_) / n_ticks_;

  const unsigned char solid_block = 219;
  const unsigned char solid_square = 254;
  for (int i = vcl_floor(progress0); i < vcl_floor(progress1); ++i)
    vcl_cout << solid_square;
  
  if (progress_ == n_ticks_)
    vcl_cout << vcl_endl;
}

#endif __NCM_TIMEBAR_H__
