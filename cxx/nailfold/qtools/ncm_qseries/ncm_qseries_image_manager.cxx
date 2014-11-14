#include "ncm_qseries_image_manager.h"

#include <vcl_iostream.h>
#include <vcl_algorithm.h>

#include <vul/vul_file.h>
#include <vul/vul_file_iterator.h>
#include <vul/vul_expand_path.h>
#include <vul/vul_string.h>

//: Define static class members here
// ncm_camera::member_ = 0;

//: Define static const class members here
// ncm_camera::const_member_;

ncm_marked_pair::ncm_marked_pair(vcl_vector< vcl_string > im_filenames, vcl_string grade_filename, bool is_graded)
: im_filenames_(im_filenames),
	grade_filename_(grade_filename),
  is_graded_(is_graded)
{
}

//: Define class member functions

//
//: Constructor
//  Generate a list of all images in the given directory, and a second list of
//  all images that do not have a corresponding markup file.
ncm_qseries_image_manager::ncm_qseries_image_manager(vcl_string root /* = "." */,
                                   vcl_string image_dir /* = "Inbox" */,
                                   vcl_string grades_dir /* = "Outbox" */)
: site_name_("unknown"),
  root_(vul_expand_path(root)),
  image_dir_(vul_expand_path(image_dir)),
  grades_dir_(vul_expand_path(grades_dir)),
  unmarked_only_(false)
{
  // You can't get return values from this but it's worth doing by default 
  // anyway - you can always run it again later
  read_site_name();
}

//
//: Read site name from <root>/site.txt and store it in site_name_
//  Return error code
//    0   Success
//    1   File not found
//    2   File empty
//    4   Unrecognized site name
int ncm_qseries_image_manager::read_site_name()
{
  // Look for a file called 'site.txt' in the root folder and read the name
  // of the site from it.
  vcl_string filename = root_ + "/site.txt";
  
  vcl_ifstream fp(filename.c_str());

  if (!fp)
    return 1;

  if (fp.eof())
    return 2;

  fp >> site_name_;

	/*char *str1 = new char[255];
	if (!fp.eof())
	{
		fp.getline(str1, 255);
		fp.getline(str1, 255);
		vcl_cout<<"Image dir read from file: " << str1 << vcl_endl;

		vcl_string im_dir(str1);
		set_image_dir(im_dir);
	}*/

  fp.close();

  // Otherwise, there's a different problem
  if (site_name() == "unknown")
    return 4;

  // Otherwise we're done
  return 0;
}

//
//: Return either a valid site name or "unknown"
vcl_string ncm_qseries_image_manager::site_name() const
{
  if (site_name_ == "ptresadern" ||
      site_name_ == "mberks" ||
      site_name_ == "ctaylor" ||
      site_name_ == "amurray" ||
      site_name_ == "gdinsdale" ||
      site_name_ == "aherrick" ||
      site_name_ == "tmoore" ||
      site_name_ == "croberts" ||

      site_name_ == "asulli" ||
      site_name_ == "mcutolo" ||
      site_name_ == "fravera" ||
      site_name_ == "vsmith" ||
      site_name_ == "rhesselstrand" ||
      site_name_ == "mwildt" ||
      site_name_ == "manderson" ||
      site_name_ == "ppyrkotsch" ||
      site_name_ == "jallen" ||
      site_name_ == "khowell" ||
      site_name_ == "hhofstee" )
    return site_name_;
  else
    return "unknown";
}

vcl_string ncm_qseries_image_manager::local_imagelist_filename() const
{
  return root_ + "/" + site_name() + "_images.txt";
}

void ncm_qseries_image_manager::set_root(vcl_string root)
{
  // Get rid of trailing slashes etc
  root_ = vul_expand_path(root);
}

void ncm_qseries_image_manager::set_image_dir(vcl_string image_dir)
{
  image_dir_ = vul_expand_path(image_dir);
}

void ncm_qseries_image_manager::set_grades_dir(vcl_string grades_dir)
{
  grades_dir_ = vul_expand_path(grades_dir);
}

int ncm_qseries_image_manager::read_image_list()
{
	//Start off with an empty image pairs list

	// Parse list of images and try to pull each from the server
  vcl_ifstream fp(local_imagelist_filename().c_str());
  if (fp == NULL)
  {
    // Return error code -1 if local image list could not be read
    return -1;
  }

  int n_failed = 0;
  int n_succeeded = 0;
  while (!fp.eof())
  {

    vcl_string filename1, filename2;
    fp >> filename1;
		fp >> filename2;
    
		vcl_vector< vcl_string> im_filenames;
		im_filenames.push_back(root_ + "/" + image_dir_ + "/" + filename1);
		im_filenames.push_back(root_ + "/" + image_dir_ + "/" + filename2);    

		vcl_string grade_filename = root_ + "/" + grades_dir_ + "/" 
			+ vul_file::strip_extension(filename1) + "_grades.txt";

		bool ims_exist = (vul_file::exists(im_filenames[0]) && vul_file::exists(im_filenames[1]));
		
    if (ims_exist) 
		{
			//Image exists, increment success count
        ++n_succeeded;

			//Check if a grade exists for the image and push onto images list
			bool is_graded = (vul_file::exists(grade_filename));
			//if (!is_graded)
			//	grade_filename = "";
			pairs_.push_back(ncm_marked_pair(im_filenames, grade_filename, is_graded));
		}
    else
      ++n_failed;
  }
  fp.close();

  // Set current image to first in the list (if applicable).
  first();

	return n_failed;
}


void ncm_qseries_image_manager::use_unmarked_only(bool unmarked_only /* = true */)
{
  unmarked_only_ = unmarked_only;
}

//
//: Return true if only showing unmarked images
bool ncm_qseries_image_manager::is_using_unmarked_only() const
{
  return unmarked_only_;
}

//
//: Count number of images in the folder, whether marked or unmarked
unsigned ncm_qseries_image_manager::n_images() const
{
  return pairs_.size();
}

//
//: Count number of unmarked images in the folder
unsigned ncm_qseries_image_manager::n_unmarked() const
{
  unsigned n_unmarked = 0;

  vcl_vector<ncm_marked_pair>::const_iterator it;
  for (it = pairs_.begin(); it != pairs_.end(); ++it)
  {
    if ((*it).is_graded_ == false)
      ++n_unmarked;
  }
  return n_unmarked;
}

//
//: Count number of valid images in the folder
unsigned ncm_qseries_image_manager::n_valid() const
{
  unsigned n_valid = 0;

  vcl_vector<ncm_marked_pair>::const_iterator it;
  for (it = pairs_.begin(); it != pairs_.end(); ++it)
  {
    if (is_valid(*it))
      ++n_valid;
  }
  return n_valid;
}

unsigned ncm_qseries_image_manager::count() const
{
  return n_valid();
}

//
//: Return true if there are no valid images in the folder
bool ncm_qseries_image_manager::is_empty() const
{
  return (count() == 0);
}

//
//: Return a list of names of valid images
vcl_vector< vcl_vector<vcl_string> > ncm_qseries_image_manager::valid_images()
{
  vcl_vector< vcl_vector<vcl_string> > images;
  
  vcl_vector<ncm_marked_pair>::iterator it = pairs_.begin();
  for (it = pairs_.begin(); it != pairs_.end(); ++it)
  {
    if (is_valid(*it))
      images.push_back((*it).im_filenames_);
  }

  return images;
}

//
//: Return the filename of the current image
vcl_vector< vcl_string > ncm_qseries_image_manager::current_image()
{
  if (pair_iterator_ != pairs_.end())
    return (*pair_iterator_).im_filenames_;
	else {
		vcl_vector< vcl_string > empty_vec;
    return empty_vec;
	}
}

//
//: Return the filename of the grade associated with the current image
vcl_string ncm_qseries_image_manager::current_grade()
{
  if (pair_iterator_ != pairs_.end())
    return (*pair_iterator_).grade_filename_;
	else {
		return "";
	}
}

//
//: Use given filename as current image. Return false if filename not in list.
bool ncm_qseries_image_manager::set_current(vcl_string filename)
{
  vcl_vector<ncm_marked_pair>::iterator it;
  for (it = pairs_.begin(); it != pairs_.end(); ++it)
  {
    // Update pair_iterator_ and return true if match found.
    if ((*it).im_filenames_[0] == filename)
    {
      pair_iterator_ = it;
      return true;
    }
  }

  // If not found, leave pair_iterator_ unchanged and return false.
  return false;
}

//
//: Tag the current image as 'marked'
void ncm_qseries_image_manager::set_current_marked()
{
  if (pair_iterator_ != pairs_.end())
    (*pair_iterator_).is_graded_ = true;
}

//
//: Return true if current image is marked
bool ncm_qseries_image_manager::is_current_marked() const
{
  return ((*pair_iterator_).is_graded_);
}

//
//: Return true if this is the first valid image in the folder
bool ncm_qseries_image_manager::is_first()
{
  vcl_vector<ncm_marked_pair>::iterator prev_image = 
      find_previous(pair_iterator_);
  return (prev_image == pairs_.end());
}

//
//: Return true if this is the last valid image in the folder
bool ncm_qseries_image_manager::is_last()
{
  vcl_vector<ncm_marked_pair>::iterator next_image = 
      find_next(pair_iterator_);
  return (next_image == pairs_.end());
}

//
//: Skip to first valid image in the folder
//  Return true if a valid entry has been found
//  pair_iterator_ is updated regardless
bool ncm_qseries_image_manager::first() 
{
  pair_iterator_ = find_first();
  return (pair_iterator_ != pairs_.end());
}
//
//: Skip to last valid image in the folder
//  Return true if a valid entry has been found
//  pair_iterator_ is updated regardless
bool ncm_qseries_image_manager::last() 
{
  pair_iterator_ = find_last();
  return (pair_iterator_ != pairs_.end());
}

//
//: Skip to previous valid image in the folder
//  Return true if a valid entry has been found
//  pair_iterator_ is only updated if a valid entry is found, otherwise it
//  remains pointing to the original image
bool ncm_qseries_image_manager::previous() 
{
  if (is_first())
    return false;
  else
  {
    pair_iterator_ = find_previous(pair_iterator_);
    return true;
  }
}

//
//: Skip to the next valid image in the folder
//  Return true if a valid entry has been found
//  pair_iterator_ is only updated if a valid entry is found, otherwise it
//  remains pointing to the original image
bool ncm_qseries_image_manager::next() 
{
  if (is_last())
    return false;
  else
  {
    pair_iterator_ = find_next(pair_iterator_);
    return true;
  }
}

//
//  Private methods
//

bool ncm_qseries_image_manager::copy_file(const vcl_string& src,
                                const vcl_string& dest) const
{
  vcl_string src_quoted = "\"" + dos_filename(src) + "\"";
  vcl_string dest_quoted = "\"" + dos_filename(dest) + "\"";
  vcl_string command_string = "copy /b " + src_quoted + " " + dest_quoted;
  bool success = (system(command_string.c_str()) == 0);

  return success;
}

bool ncm_qseries_image_manager::move_file(const vcl_string& src,
                                const vcl_string& dest) const
{
  vcl_string src_quoted = "\"" + dos_filename(src) + "\"";
  vcl_string dest_quoted = "\"" + dos_filename(dest) + "\"";
  vcl_string command_string = "move " + src_quoted + " " + dest_quoted;
  bool success = (system(command_string.c_str()) == 0);

  return success;
}


vcl_string ncm_qseries_image_manager::dos_filename(const vcl_string& filename) const
{
  vcl_string dos_filename = filename;
  vul_string_replace(dos_filename, "/", "\\");
  return dos_filename;
}

//
//: Check validity of a given marked iterator
bool ncm_qseries_image_manager::is_valid(const ncm_marked_pair& marked_image) const
{
  // If we are only interested in unmarked images then check for this, too
  if (unmarked_only_)
    return (marked_image.is_graded_ == false);
  else
    return true;
}

//
//: Find next valid iterator after 'initial'
//  Return pairs_.end() if none are valid
vcl_vector<ncm_marked_pair>::iterator ncm_qseries_image_manager::find_next(
    const vcl_vector<ncm_marked_pair>::iterator& initial)
{
  vcl_vector<ncm_marked_pair>::iterator it = initial;

  ++it; // skip current element in vector
  while (it < pairs_.end())
  {
    if (is_valid(*it))
      return it;
    else
      ++it;
  }

  return pairs_.end();
}

//
//: Find previous valid iterator
//  Return pairs_.end() if none are valid
vcl_vector<ncm_marked_pair>::iterator ncm_qseries_image_manager::find_previous(
    const vcl_vector<ncm_marked_pair>::iterator& initial)
{
  vcl_vector<ncm_marked_pair>::reverse_iterator rit(initial);

  // No need to decrement here as we are using a reverse-iterator which
  // automatically looks 'backwards'.
  while (rit < pairs_.rend())
  {
    if (is_valid(*rit))
      // But when converting from a reverse- to forward-iterator, we need to
      // offset by 1 
      return vcl_vector<ncm_marked_pair>::iterator(rit.base()-1);
    else
      ++rit;
  }

  return pairs_.end();
}

//
//: Find first valid iterator
//  Return pairs_.end() if none are valid
vcl_vector<ncm_marked_pair>::iterator ncm_qseries_image_manager::find_first()
{
  if (is_empty())
    return pairs_.end();
  else if (is_valid(*pairs_.begin()))
    return pairs_.begin();
  else
    return find_next(pairs_.begin());
}

//
//: Find last valid iterator
//  Return pairs_.end() if none are valid
vcl_vector<ncm_marked_pair>::iterator ncm_qseries_image_manager::find_last()
{
  if (is_empty())
    return pairs_.end();
  else
    return find_previous(pairs_.end());
}
