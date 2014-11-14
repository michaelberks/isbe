#include "ncm_qmarkup_server.h"

#include <vcl_iostream.h>
#include <vcl_algorithm.h>

#include <vul/vul_file.h>
#include <vul/vul_file_iterator.h>
#include <vul/vul_expand_path.h>
#include <vul/vul_string.h>

#include <nailfold/ncm_server_handler.h>

//: Define static class members here
// ncm_camera::member_ = 0;

//: Define static const class members here
// ncm_camera::const_member_;

ncm_marked_image::ncm_marked_image(vcl_string filename, bool is_marked)
: filename_(filename),
  is_marked_(is_marked)
{
}

//: Define class member functions

//
//: Constructor
//  Generate a list of all images in the given directory, and a second list of
//  all images that do not have a corresponding markup file.
ncm_qmarkup_sever::ncm_qmarkup_sever(vcl_string root /* = "." */,
                                   vcl_string inbox /* = "Inbox" */,
                                   vcl_string outbox /* = "Outbox" */,
                                   vcl_string sent /* = "Sent" */)
: site_name_("unknown"),
  root_(vul_expand_path(root)),
  inbox_(vul_expand_path(inbox)),
  outbox_(vul_expand_path(outbox)),
  sent_(vul_expand_path(sent)),
  unmarked_only_(true)
{
  parse_inbox();

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
int ncm_qmarkup_sever::read_site_name()
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
  fp.close();

  // Otherwise, there's a different problem
  if (site_name() == "unknown")
    return 4;

  // Otherwise we're done
  return 0;
}

//
//: Return either a valid site name or "unknown"
vcl_string ncm_qmarkup_sever::site_name() const
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

vcl_string ncm_qmarkup_sever::remote_imagelist_filename() const
{
  return site_name() + "_images.txt";
}
vcl_string ncm_qmarkup_sever::local_imagelist_filename() const
{
  return root_ + "/" + remote_imagelist_filename();
}

void ncm_qmarkup_sever::set_root(vcl_string root)
{
  // Get rid of trailing slashes etc
  root_ = vul_expand_path(root);
  parse_inbox();
}

void ncm_qmarkup_sever::set_inbox(vcl_string inbox)
{
  inbox_ = inbox;
  parse_inbox();
}

void ncm_qmarkup_sever::set_outbox(vcl_string outbox)
{
  outbox_ = vul_expand_path(outbox);
}

void ncm_qmarkup_sever::set_sent(vcl_string sent)
{
  sent_ = vul_expand_path(sent);
}

void ncm_qmarkup_sever::parse_inbox()
{
  vcl_vector<vcl_string> extensions;
  extensions.push_back("jpg");
  extensions.push_back("png");
  extensions.push_back("bmp");

  images_.clear();

  for (unsigned i = 0; i < extensions.size(); ++i)
  {
    vcl_string search_string = root_ + "/" + inbox_ + "/*." + extensions[i];
    for (vul_file_iterator fn = search_string; fn; ++fn) 
    {
      const vcl_string markup_filename = 
          vul_file::strip_extension(fn()) + "_markup.txt";

      const bool is_marked = (vul_file::exists(markup_filename));
      images_.push_back(ncm_marked_image(fn(), is_marked));
    }
  }

  // Set current image to first in the list (if applicable).
  first();
}

//
//: Retrieve any missing image files from the remote server to the Inbox.
//    n_images = -1 => download image list, grab all images and parse inbox
//    n_images = 0  => download the image list only
//    n_images > 0  => download n_images images only (using the existing 
//                     image list)
int ncm_qmarkup_sever::pull_imagelist_from(const ncm_server_handler& server_handler)
{
  // Get list of images from the remote server
  const bool success = server_handler.pull(remote_imagelist_filename(), 
                                           local_imagelist_filename());

  // Return error code -1 if failed to grab the image list for the site
  if (!success)
    return -1;
  
  return 0;
}

//
//: Retrieve any missing image files from the remote server to the Inbox.
//    n_images = -1 => grab all images
//    n_images >= 0 => download n_images images only (using the existing 
//                     image list)
int ncm_qmarkup_sever::pull_images_from(const ncm_server_handler& server_handler,
                                      int n_images /* = -1 */)
{
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
    if ((n_images >= 0) && (n_succeeded == n_images))
      break;

    vcl_string filename;
    fp >> filename;

    vcl_string dest = root_ + "/" + inbox_ + "/" + filename;

    // Ignore files that already exist
    if (!vul_file::exists(dest))
    {
      vcl_string src = "images/" + filename;
      const bool image_success = server_handler.pull(src, dest);

      if (image_success)
        ++n_succeeded;
      else
        ++n_failed;
    }
  }
  fp.close();

  return n_failed;
}

//
//: Retrieve any existing markup files from the remote server to the Sent 
//  folder, then transfer most recent markup files to the Inbox.
int ncm_qmarkup_sever::pull_markup_from(const ncm_server_handler& server_handler)
{
  vcl_string src = "markup/" + site_name() + "/*_markup.txt";
  vcl_string dest = root_ + "/" + sent_ + "/";

  // Additional options (-unsafe) are needed here to use server-side wildcards
  int success = server_handler.pull(src, dest, "mµî³¡æ¥");

  if (!success)
    return 1;

  // Get a list of all completed markups (i.e. everything within Sent).
  // (We could also add the files from Outbox here to be doubly sure.)
  vcl_vector<vcl_string> markup_filenames;
  vcl_string search_string = root_ + "/" + sent_ + "/*_markup.txt";
  for (vul_file_iterator fn = search_string; fn; ++fn) 
    markup_filenames.push_back(fn());

  // Sort the filenames so that they are in date order (as files in Sent are
  // timestamped)
  vcl_sort(markup_filenames.begin(), markup_filenames.end());

  // Now copy files to Inbox in reverse order, ignoring anything that already 
  // exists (such that the newest file is copied first and older files are
  // ignored).
  vcl_vector<vcl_string>::reverse_iterator rit;
  for (rit = markup_filenames.rbegin(); rit != markup_filenames.rend(); ++rit)
  {
    unsigned start_char = (*rit).find('#');
    const vcl_string filename = (*rit).substr(start_char+1);
    const vcl_string dest_markup = root_ + "/" + inbox_ + "/" + filename;

    if (!vul_file::exists(dest_markup))
      copy_file(*rit, dest_markup);
  }

  return 0;
}

//
//: Copy named file from Inbox to Outbox
bool ncm_qmarkup_sever::copy_to_outbox(vcl_string src,
                                     vcl_string dest /* = "" */)
{
  // If no destination filename given then use source filename
  if (dest == "")
    dest = src;

  // Source and destination file names
  // Surrounded by quotes to accommodate folder and file names with spaces
  src = root_ + "/" + inbox_ + "/" + vul_file::strip_directory(src);
  dest = root_ + "/" + outbox_ + "/" + vul_file::strip_directory(dest);

  return copy_file(src, dest);
}

//
//: Copy any markup files in Outbox to remote server and move the local copies
//  to the Sent folder.
int ncm_qmarkup_sever::push_outbox_to(const ncm_server_handler& server_handler)
{
  vcl_string search_string = root_ + "/" + outbox_ + "/*_markup.txt";
  int n_failed = 0;
  for (vul_file_iterator fn = search_string; fn; ++fn)
  {
    vcl_string src = fn();
    vcl_string dest = "markup/" + site_name() + "/" + fn.filename();
    bool push_success = server_handler.push(src, dest);

    if (push_success)
    {
      // Source and destination file names
      // Surrounded by quotes to accommodate folder and file names with spaces
      vcl_string src = fn();
      vcl_string dest = root_ + "/" + sent_ + "/" + fn.filename();
      move_file(src, dest);
    }
    else
    {
      ++n_failed;
    }
  }

  return n_failed;
}


void ncm_qmarkup_sever::use_unmarked_only(bool unmarked_only /* = true */)
{
  unmarked_only_ = unmarked_only;
}

//
//: Return true if only showing unmarked images
bool ncm_qmarkup_sever::is_using_unmarked_only() const
{
  return unmarked_only_;
}

//
//: Count number of images in the folder, whether marked or unmarked
unsigned ncm_qmarkup_sever::n_images() const
{
  return images_.size();
}

//
//: Count number of unmarked images in the folder
unsigned ncm_qmarkup_sever::n_unmarked() const
{
  unsigned n_unmarked = 0;

  vcl_vector<ncm_marked_image>::const_iterator it;
  for (it = images_.begin(); it != images_.end(); ++it)
  {
    if ((*it).is_marked_ == false)
      ++n_unmarked;
  }
  return n_unmarked;
}

//
//: Count number of valid images in the folder
unsigned ncm_qmarkup_sever::n_valid() const
{
  unsigned n_valid = 0;

  vcl_vector<ncm_marked_image>::const_iterator it;
  for (it = images_.begin(); it != images_.end(); ++it)
  {
    if (is_valid(*it))
      ++n_valid;
  }
  return n_valid;
}

unsigned ncm_qmarkup_sever::count() const
{
  return n_valid();
}

//
//: Return true if there are no valid images in the folder
bool ncm_qmarkup_sever::is_empty() const
{
  return (count() == 0);
}

//
//: Return a list of names of valid images
vcl_vector<vcl_string> ncm_qmarkup_sever::valid_images()
{
  vcl_vector<vcl_string> images;
  
  vcl_vector<ncm_marked_image>::iterator it = images_.begin();
  for (it = images_.begin(); it != images_.end(); ++it)
  {
    if (is_valid(*it))
      images.push_back((*it).filename_);
  }

  return images;
}

//
//: Return the filename of the current image
vcl_string ncm_qmarkup_sever::current_image()
{
  if (image_iterator_ != images_.end())
    return (*image_iterator_).filename_;
  else
    return "";
}

//
//: Use given filename as current image. Return false if filename not in list.
bool ncm_qmarkup_sever::set_current(vcl_string filename)
{
  vcl_vector<ncm_marked_image>::iterator it;
  for (it = images_.begin(); it != images_.end(); ++it)
  {
    // Update image_iterator_ and return true if match found.
    if ((*it).filename_ == filename)
    {
      image_iterator_ = it;
      return true;
    }
  }

  // If not found, leave image_iterator_ unchanged and return false.
  return false;
}

//
//: Tag the current image as 'marked'
void ncm_qmarkup_sever::set_current_marked()
{
  if (image_iterator_ != images_.end())
    (*image_iterator_).is_marked_ = true;
}

//
//: Return true if current image is marked
bool ncm_qmarkup_sever::is_current_marked() const
{
  return ((*image_iterator_).is_marked_);
}

//
//: Return true if this is the first valid image in the folder
bool ncm_qmarkup_sever::is_first()
{
  vcl_vector<ncm_marked_image>::iterator prev_image = 
      find_previous(image_iterator_);
  return (prev_image == images_.end());
}

//
//: Return true if this is the last valid image in the folder
bool ncm_qmarkup_sever::is_last()
{
  vcl_vector<ncm_marked_image>::iterator next_image = 
      find_next(image_iterator_);
  return (next_image == images_.end());
}

//
//: Skip to first valid image in the folder
//  Return true if a valid entry has been found
//  image_iterator_ is updated regardless
bool ncm_qmarkup_sever::first() 
{
  image_iterator_ = find_first();
  return (image_iterator_ != images_.end());
}
//
//: Skip to last valid image in the folder
//  Return true if a valid entry has been found
//  image_iterator_ is updated regardless
bool ncm_qmarkup_sever::last() 
{
  image_iterator_ = find_last();
  return (image_iterator_ != images_.end());
}

//
//: Skip to previous valid image in the folder
//  Return true if a valid entry has been found
//  image_iterator_ is only updated if a valid entry is found, otherwise it
//  remains pointing to the original image
bool ncm_qmarkup_sever::previous() 
{
  if (is_first())
    return false;
  else
  {
    image_iterator_ = find_previous(image_iterator_);
    return true;
  }
}

//
//: Skip to the next valid image in the folder
//  Return true if a valid entry has been found
//  image_iterator_ is only updated if a valid entry is found, otherwise it
//  remains pointing to the original image
bool ncm_qmarkup_sever::next() 
{
  if (is_last())
    return false;
  else
  {
    image_iterator_ = find_next(image_iterator_);
    return true;
  }
}

//
//  Private methods
//

bool ncm_qmarkup_sever::copy_file(const vcl_string& src,
                                const vcl_string& dest) const
{
  vcl_string src_quoted = "\"" + dos_filename(src) + "\"";
  vcl_string dest_quoted = "\"" + dos_filename(dest) + "\"";
  vcl_string command_string = "copy /b " + src_quoted + " " + dest_quoted;
  bool success = (system(command_string.c_str()) == 0);

  return success;
}

bool ncm_qmarkup_sever::move_file(const vcl_string& src,
                                const vcl_string& dest) const
{
  vcl_string src_quoted = "\"" + dos_filename(src) + "\"";
  vcl_string dest_quoted = "\"" + dos_filename(dest) + "\"";
  vcl_string command_string = "move " + src_quoted + " " + dest_quoted;
  bool success = (system(command_string.c_str()) == 0);

  return success;
}


vcl_string ncm_qmarkup_sever::dos_filename(const vcl_string& filename) const
{
  vcl_string dos_filename = filename;
  vul_string_replace(dos_filename, "/", "\\");
  return dos_filename;
}

//
//: Check validity of a given marked iterator
bool ncm_qmarkup_sever::is_valid(const ncm_marked_image& marked_image) const
{
  // If we are only interested in unmarked images then check for this, too
  if (unmarked_only_)
    return (marked_image.is_marked_ == false);
  else
    return true;
}

//
//: Find next valid iterator after 'initial'
//  Return images_.end() if none are valid
vcl_vector<ncm_marked_image>::iterator ncm_qmarkup_sever::find_next(
    const vcl_vector<ncm_marked_image>::iterator& initial)
{
  vcl_vector<ncm_marked_image>::iterator it = initial;

  ++it; // skip current element in vector
  while (it < images_.end())
  {
    if (is_valid(*it))
      return it;
    else
      ++it;
  }

  return images_.end();
}

//
//: Find previous valid iterator
//  Return images_.end() if none are valid
vcl_vector<ncm_marked_image>::iterator ncm_qmarkup_sever::find_previous(
    const vcl_vector<ncm_marked_image>::iterator& initial)
{
  vcl_vector<ncm_marked_image>::reverse_iterator rit(initial);

  // No need to decrement here as we are using a reverse-iterator which
  // automatically looks 'backwards'.
  while (rit < images_.rend())
  {
    if (is_valid(*rit))
      // But when converting from a reverse- to forward-iterator, we need to
      // offset by 1 
      return vcl_vector<ncm_marked_image>::iterator(rit.base()-1);
    else
      ++rit;
  }

  return images_.end();
}

//
//: Find first valid iterator
//  Return images_.end() if none are valid
vcl_vector<ncm_marked_image>::iterator ncm_qmarkup_sever::find_first()
{
  if (is_empty())
    return images_.end();
  else if (is_valid(*images_.begin()))
    return images_.begin();
  else
    return find_next(images_.begin());
}

//
//: Find last valid iterator
//  Return images_.end() if none are valid
vcl_vector<ncm_marked_image>::iterator ncm_qmarkup_sever::find_last()
{
  if (is_empty())
    return images_.end();
  else
    return find_previous(images_.end());
}
