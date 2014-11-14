#include "ncm_image_grade.h"

#include <nailfold/ncm_annotation.h>

//: Default constructor
ncm_image_grade::ncm_image_grade(ncm_annotation* parent_annotation /* = NULL */)
: image_grade_(GradeUndefined),
  parent_annotation_(parent_annotation)
{
}

//: Destructor
ncm_image_grade::~ncm_image_grade()
{
}

void ncm_image_grade::set_value(ImageGrade image_grade)
{
  if (image_grade_ == image_grade)
    return;

  image_grade_ = image_grade;

  if (parent_annotation_ != NULL)
    parent_annotation_->set_modified();
}

void ncm_image_grade::set_value(int image_grade)
{
  set_value(static_cast<ImageGrade>(image_grade));
}

ncm_image_grade::ImageGrade ncm_image_grade::value() const
{
  return image_grade_;
}

//
//: Image grade as a string
void ncm_image_grade::set_from_string(vcl_string image_grade)
{
  // Should probably reimplement this as a map

  if (image_grade == "Undefined")
    set_value(GradeUndefined);
  else if (image_grade == "Normal")
    set_value(GradeNormal);
  else if (image_grade == "Early")
    set_value(GradeEarly);
  else if (image_grade == "Active")
    set_value(GradeActive);
  else if (image_grade == "Late")
    set_value(GradeLate);
  else if (image_grade == "Ungradeable_Condition")
    set_value(GradeExtreme);
  else if (image_grade == "Ungradeable_Quality")
    set_value(GradePoorQuality);
  else if (image_grade == "Non-specific")
    set_value(GradeNonspecific);
	else if (image_grade == "auto")
    set_value(GradeUndefined);
  else
    assert(false);
}

vcl_string ncm_image_grade::as_string() const
{
  switch (image_grade_)
  {
    case GradeUndefined: return "Undefined";
    case GradeNormal: return "Normal";
    case GradeEarly: return "Early";
    case GradeActive: return "Active";
    case GradeLate: return "Late";
    case GradeExtreme: return "Ungradeable_Condition";
    case GradePoorQuality: return "Ungradeable_Quality";
    case GradeNonspecific: return "Non-specific";
    default: assert(false);
  }
  return "";
}
