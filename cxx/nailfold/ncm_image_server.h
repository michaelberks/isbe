#ifndef ncm_file_server_h_
#define ncm_file_server_h_

#include <vcl_string.h>
#include <vcl_vector.h>

class ncm_server_handler;
class ncm_marked_image;

class ncm_file_server
{
// INTERFACE

public:
  // No member variables here, please

  //: Default ctor
  ncm_file_server(vcl_string root = ".",
                  vcl_string inbox = "Inbox",
                  vcl_string outbox = "Outbox",
                  vcl_string sent = "Sent");

  // The 'Big Three'

  ////: Copy ctor
  //ncm_camera(ncm_camera const& rhs);
  ////: dtor
  //~ncm_camera();
  ////: Assignment operator
  //ncm_camera& operator=(ncm_camera const& rhs); // assignment

  int read_site_name();
  vcl_string site_name() const;

  vcl_string remote_imagelist_filename() const;
  vcl_string local_imagelist_filename() const;

  void set_root(vcl_string root);
	void set_user(vcl_string user);
  void set_image_dir(vcl_string image_dir);
  void set_markup_dir(vcl_string markup_dir);
  void set_inbox(vcl_string inbox);
  void set_outbox(vcl_string outbox);
  void set_sent(vcl_string sent);

  void parse_inbox();
	void parse_images();
  vcl_vector<vcl_string> valid_images();

  unsigned n_images() const;
  unsigned n_unmarked() const;
  unsigned n_valid() const;
  unsigned count() const;
  bool is_empty() const;

  void use_unmarked_only(bool unmarked_only = true);
  bool is_using_unmarked_only() const;

  vcl_string current_image();
  void set_current_marked();
  bool is_current_marked() const;
  bool set_current(vcl_string filename);

  int pull_imagelist_from(const ncm_server_handler& server_handler);
  int pull_images_from(const ncm_server_handler& server_handler,
                       int n_images = -1);

  int pull_markup_from(const ncm_server_handler& server_handler);

  int push_outbox_to(const ncm_server_handler& server_handler);
	int push_outbox_to(const vcl_string local_path);

  bool copy_to_outbox(vcl_string src,
                      vcl_string dest = "");

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

  bool is_valid(const ncm_marked_image& marked_image) const;
  
  // Get iterators for first, last, previous and next item
  vcl_vector<ncm_marked_image>::iterator find_previous(
      const vcl_vector<ncm_marked_image>::iterator& initial);
  vcl_vector<ncm_marked_image>::iterator find_next(
      const vcl_vector<ncm_marked_image>::iterator& initial);
  vcl_vector<ncm_marked_image>::iterator find_first();
  vcl_vector<ncm_marked_image>::iterator find_last();

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
	vcl_string user_;
  vcl_string image_dir_;
  vcl_string markup_dir_;
  vcl_string inbox_;
  vcl_string outbox_;
  vcl_string sent_;

  vcl_vector<ncm_marked_image> images_;
  vcl_vector<ncm_marked_image>::iterator image_iterator_;
};


//: Helper class that stores the image filename with the marked status
class ncm_marked_image 
{
  friend ncm_file_server;

public:
  ncm_marked_image(vcl_string filename, bool is_marked);

private:
  vcl_string filename_;
  bool is_marked_;
};

#endif // ncm_file_server_h_