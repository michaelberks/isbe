#include "nailfold/qtools/ncm_qseries/ncm_qseries_grade.h"

#include <vcl_ctime.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>
#include <vcl_sstream.h>
#include <vcl_iomanip.h>
#include <vcl_sstream.h>
#include <vcl_iomanip.h>

#include <vul/vul_file.h>
#include <vul/vul_string.h>

#include <vsl/vsl_binary_io.h>
//#include <vsl/vsl_binary_loader.h>
#include <vsl/vsl_indent.h>

#include <mbl/mbl_parse_block.h>
#include <mbl/mbl_read_props.h>
#include <mbl/mbl_parse_keyword_list.h>
#include <mbl/mbl_index_sort.h>

#include <nailfold/ncm_encrypt.h>
#include <nailfold/ncm_vessel.h>
#include <nailfold/ncm_vessel_properties.h>
#include <nailfold/ncm_haemorrhage.h>

//: Define static class members here
// type ncm_qseries_grade::member_ = 0;

//: Define static const class members here
// type ncm_qseries_grade::const_member_ = 0;

//
// Public methods
//

//
//: Default constructor
ncm_qseries_grade::ncm_qseries_grade()
: snapshot_on_modify_(true)
{
	filename_ = "";
  clear();
}

//
//: Destructor
ncm_qseries_grade::~ncm_qseries_grade()
{
  set_snapshots(false);
  clear();
}

//
//: Clear all data
void ncm_qseries_grade::clear()
{
  observer_name_ = "-";
  datetime_ = static_cast<vcl_time_t>(0);
  time_taken_ = 0.0;

	//Initialise grades and reasons
	is_modified_ = false;
	is_graded_ = false;
	progression_level_ = None;
	progress_reasons_.resize( NumReasons + 1 );

	//Initialise to false just in case...
	for ( vcl_vector<bool>::iterator itr = progress_reasons_.begin(), end = progress_reasons_.end(); itr != end; ++itr )
		 *itr = false;
  
}

//
//: Markup filename
void ncm_qseries_grade::set_filename(const vcl_string& filename)
{
  filename_ = filename;
}
vcl_string ncm_qseries_grade::filename() const
{
  return filename_;
}

//
//: Generate a tagged copy of the filename in the format:
//  YYMMDD-hhmmss_<observer_name>#<filename>
vcl_string ncm_qseries_grade::tagged_filename() const
{
  // Create timestamp
  time_t rawtime = time(NULL);
  struct tm* ptm = gmtime(&rawtime);
  vcl_stringstream ss;
  ss << ptm->tm_year+1900 
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_mon+1 // months start at 0
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_mday
     << "-"
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_hour 
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_min 
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_sec;
  vcl_string timestamp = ss.str();

  return timestamp + "_" + 
         observer_name_ + "#" + 
         vul_file::strip_directory(filename_);
}

bool ncm_qseries_grade::save()
{
  if (filename_.empty())
    return false;

  // else
  stop_timing();
  vcl_ofstream ofs(filename_.c_str());
  t_write(ofs);
  ofs.close();

  // No longer modified since last save
  is_modified_ = false;

  return true;
}

void ncm_qseries_grade::set_snapshots(const bool snapshot_on_modify /* = true */)
{
  snapshot_on_modify_ = snapshot_on_modify;
}
bool ncm_qseries_grade::snapshot() const
{
  // else
  vcl_ofstream ofs(snapshot_filename().c_str());
  t_write(ofs);
  ofs.close();
  return true;
}

vcl_string ncm_qseries_grade::snapshot_filename() const
{
  return filename_ + "~";
}


//
//: Set the observer's name
void ncm_qseries_grade::set_observer_name(const vcl_string& observer_name)
{
  if (observer_name_ != observer_name)
  {
    observer_name_ = observer_name;
    set_modified();
  }
}

//
//: Reset clock timer (used for storing annotation time)
void ncm_qseries_grade::start_timing()
{
  vcl_time(&datetime_);

  // Don't reset time_taken_ as we will add the time spent annotating to the
  // stored time (if any) to reflect the total time spent marking up the image

  // Don't use this as a cue to set_modified() - instead, we'll wait until the 
  // user makes an observable change to the markup.
}

//
//: Store length of time taken to annotate image
void ncm_qseries_grade::stop_timing()
{
  // Update total time spent so far marking up this image
  time_taken_ += vcl_difftime(vcl_time(NULL), datetime_);
  //set_modified();
}

//
//: Return date and time as a string
vcl_string ncm_qseries_grade::date_string() const
{
  char* date_string = vcl_asctime(vcl_gmtime(&datetime_));

  if (date_string != NULL)
    return date_string;
  else
    return "-";
}

//
//: Annotation has been modified (or force reset)
void ncm_qseries_grade::set_modified(const bool is_modified /* = true */) const
{
  // is_modified_ is mutable and therefore allowed to change
  is_modified_ = is_modified;

  if (is_modified && snapshot_on_modify_)
    snapshot();
}

void ncm_qseries_grade::set_normal(const bool is_normal)
{
  // is_modified_ is mutable and therefore allowed to change
  is_normal_ = is_normal;
	is_graded_ = true;
	set_modified(true);
	vcl_cout << "Are the images are normal? User says: " << is_normal << vcl_endl;
}

void ncm_qseries_grade::set_reason(const bool selected, const ProgressReason reason)
{
	progress_reasons_[ reason ] = selected;
	//vcl_vector<bool>::iterator itr = progress_reasons_.begin();
	//*(itr+reason) = selected;
	vcl_cout << "Reason " << reason << " selected as " << selected << vcl_endl;
	set_modified(true);
}

void ncm_qseries_grade::set_progression(const ProgressionLevel level)
{
	progression_level_ = level;
}

//
//: Whether the annotation has been modified
bool ncm_qseries_grade::is_modified() const
{
  return is_modified_;
}

 
//: Return union of missing data flags to indicate missing annotations.
bool ncm_qseries_grade::is_graded() const
{
  return is_graded_;
}

//Are the images marked as being normal
bool ncm_qseries_grade::is_normal() const
{
  return is_normal_;
}

//Have any reasons been selected
bool ncm_qseries_grade::is_reasons()
{
	for ( vcl_vector<bool>::iterator itr = progress_reasons_.begin(), end = progress_reasons_.end(); itr != end; ++itr )
	{
		if (*itr)
			return true;
	}
	return false;
}

ncm_qseries_grade::ProgressionLevel ncm_qseries_grade::progression_level()
{
	return progression_level_;
}

//Get selected reasons
vcl_vector< bool > ncm_qseries_grade::progress_reasons()
{
	return progress_reasons_;
}
bool ncm_qseries_grade::progress_reasons(ProgressReason reason)
{
	return progress_reasons_[reason];
	//vcl_vector<bool>::iterator itr = progress_reasons_.begin();
	//return *(itr+reason);
}

//
//: Version number of this class's binary IO
short ncm_qseries_grade::version() const
{ 
  return 2; // QMarkup v1.4
}

//
//: The class name
vcl_string ncm_qseries_grade::is_a() const
{
  return "ncm_qseries_grade"; 
}

//
//: Class name equivalence
bool ncm_qseries_grade::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str);
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_qseries_grade& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_qseries_grade& b)
{
  b.b_read(bfs);
}

//
//: Text write
void ncm_qseries_grade::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {" << '\n';

  vsl_indent_inc(tfs);
  tfs << vsl_indent() << "version: " << version() << '\n';

  //// version 3
  //vsl_b_write(bfs,new_variable_);

  // version 2
    // No change - just used as a marker to keep track of the qmarkup 
    // software being used.

  // version 1
  tfs << vsl_indent() << "observer: " << observer_name_ << '\n';
  tfs << vsl_indent() << "timestamp: " << date_string(); // EOL appended already
  tfs << vsl_indent() << "time_taken: " << time_taken_ << '\n';
	tfs << vsl_indent() << "ims_normal: " << is_normal_ << '\n';

	if (!is_normal())
	{
		//Write out progression level
		tfs << vsl_indent() << "progression_level: " << progression_level_ << '\n';

		if (progression_level_ != None)
		{
			//Write out reasons for progression
			int reason = 0;
			for ( vcl_vector<bool>::const_iterator itr = progress_reasons_.begin(), end = progress_reasons_.end(); itr != end; ++itr )
			{
				tfs << vsl_indent() << "r" << reason << "_selected: " << *itr << '\n';
				reason++;
			}
		}
	}

	tfs << vsl_indent() << "} // " << is_a() << '\n';
  tfs << vcl_flush;
}

//
//: Text read (from stream, without identifying token)
void ncm_qseries_grade::t_read(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  // Disable snapshots while we load in new data
  set_snapshots(false);
  clear();

  mbl_read_props_type props = mbl_read_props(tfs);

  observer_name_ =  props.get_required_property("observer");

  // Timestamp is not read here as it will be reset to reflect the last
  // time that an image was edited

  time_taken_ = vul_string_atof(props.get_required_property("time_taken"));

	//Are the images normal
	const bool ims_normal = vul_string_atof(props.get_required_property("ims_normal")) > 0;
	set_normal(ims_normal);

	if (!is_normal())
	{
		//Was any level of progression given?
		progression_level_ = static_cast<ProgressionLevel>(vul_string_atoi(props.get_required_property("progression_level")));

		if (progression_level_ != None)
		{
			//What reasons were selected?
			int reason = 0;
			for ( vcl_vector<bool>::iterator itr = progress_reasons_.begin(), end = progress_reasons_.end(); itr != end; ++itr )
			{
				vcl_stringstream ss;
				ss << "r" << reason << "_selected";
				bool reason_set = vul_string_atof(props.get_required_property(ss.str())) > 0;
				*itr = reason_set;
				reason++;
			}
		}
	}

  is_modified_ = false;

  set_snapshots(true);
}

//
//: Text read (from file, with identifying token)
void ncm_qseries_grade::t_read(const vcl_string& filename)
{
  vcl_ifstream tfs(filename.c_str());

  mbl_read_props_type props = mbl_read_props(tfs);
  vcl_string token = is_a();
  vcl_string markup_string = props.get_required_property(token);
  vcl_stringstream markup_stream(markup_string);
  t_read(markup_stream);

  tfs.close();
}

//
//: Text write by reference
vcl_ostream& operator<<(vcl_ostream& os, const ncm_qseries_grade& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text write by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_qseries_grade* b)
{
  if (b)
    return os << *b;
  else
    return os << "No " << b->is_a() << " defined.";
}


//
// Protected functions

//
//: Print summary to output stream os
void ncm_qseries_grade::print_summary(vcl_ostream& os) const
{
  os << "Observer: " << observer_name_ << '\n';
  os << "Date: " << "<date>" << '\n';
  os << "Time taken: " << time_taken_ << '\n';
  os << vcl_flush;
}

//
//: Binary write
void ncm_qseries_grade::b_write(vsl_b_ostream& bfs) const
{
  vsl_b_write(bfs,version());

  //// version 2
  //vsl_b_write(bfs,new_variable_);

  // version 1
  vsl_b_write(bfs,observer_name_);
  vsl_b_write(bfs,datetime_);
  vsl_b_write(bfs,time_taken_);
}

//
//: Binary read
void ncm_qseries_grade::b_read(vsl_b_istream& bfs)
{
  short version;
  vsl_b_read(bfs,version);
  switch (version)
  {
    //case (2): 
    //  vsl_b_read(bfs,new_variable_);
    //  // no break

    case (1):
      vsl_b_read(bfs,observer_name_);
      vsl_b_read(bfs,datetime_);
      vsl_b_read(bfs,time_taken_);

      break;

    default:
      vcl_cerr << is_a() << "::b_read() ";
      vcl_cerr << "Unexpected version number " << version << vcl_endl;
      abort();
  }
}

//
//  Private methods
//
