#ifndef ncm_qseries_image_manager_h_
#define ncm_qseries_image_manager_h_

#include <vcl_string.h>
#include <vcl_vector.h>

class ncm_server_handler;
class ncm_marked_pair;

class ncm_qseries_image_manager
{
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_qseries_image_manager(vcl_string root = ".",
                  vcl_string image_dir = "images",
                  vcl_string grades_dir = "grades");

  // The 'Big Three'

  ////: Copy ctor
  //ncm_camera(ncm_camera const& rhs);
  ////: dtor
  //~ncm_camera();
  ////: Assignment operator
  //ncm_camera& operator=(ncm_camera const& rhs); // assignment

  int read_site_name();
  vcl_string site_name() const;

  vcl_string local_imagelist_filename() const;

  void set_root(vcl_string root);
  void set_image_dir(vcl_string inbox);
  void set_grades_dir(vcl_string outbox);

  int read_image_list();
  vcl_vector< vcl_vector<vcl_string> > valid_images();

  unsigned n_images() const;
  unsigned n_unmarked() const;
  unsigned n_valid() const;
  unsigned count() const;
  bool is_empty() const;

  void use_unmarked_only(bool unmarked_only = true);
  bool is_using_unmarked_only() const;

  vcl_vector< vcl_string > current_image();
	vcl_string current_grade();
  void set_current_marked();
  bool is_current_marked() const;
  bool set_current(vcl_string filename);

  bool is_first();
  bool is_last();

  bool first();
  bool last();
  bool previous();
  bool next();


// IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

private:
  // Members and functions visible only to objects of this class

  bool is_valid(const ncm_marked_pair& marked_pair) const;
  
  // Get iterators for first, last, previous and next item
  vcl_vector<ncm_marked_pair>::iterator find_previous(
      const vcl_vector<ncm_marked_pair>::iterator& initial);
  vcl_vector<ncm_marked_pair>::iterator find_next(
      const vcl_vector<ncm_marked_pair>::iterator& initial);
  vcl_vector<ncm_marked_pair>::iterator find_first();
  vcl_vector<ncm_marked_pair>::iterator find_last();

  bool copy_file(const vcl_string& src,
                 const vcl_string& dest) const;
  
  bool move_file(const vcl_string& src,
                 const vcl_string& dest) const;

  vcl_string dos_filename(const vcl_string& filename) const;

  //: Name of the site (e.g. manchester, liverpool, salford)
  vcl_string site_name_;

  //: Whether to consider only unmarked images
  bool unmarked_only_;

  vcl_string root_;
  vcl_string image_dir_;
  vcl_string grades_dir_;

  vcl_vector<ncm_marked_pair> pairs_;
  vcl_vector<ncm_marked_pair>::iterator pair_iterator_;
};


//: Helper class that stores the image filename with the marked status
class ncm_marked_pair 
{
  friend ncm_qseries_image_manager;

public:
  ncm_marked_pair(vcl_vector< vcl_string > im_filenames, vcl_string grade_filename, bool is_marked);

private:
  vcl_vector< vcl_string > im_filenames_;
	vcl_string grade_filename_;
  bool is_graded_;
};

#endif // ncm_qseries_image_manager_h_