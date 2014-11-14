#ifndef __ncm_qvessel_properties_h__
#define __ncm_qvessel_properties_h__

//:
// \file
// \brief Extension of ncm_vessel_properties to accept signals from Qt objects
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <QObject>

#include <nailfold/ncm_vessel_properties.h>

//  To enable signals and slots, the object must be a derived class of QObject.
//  To do this I've used multiple inheritance, though this is known to be
//  not for the fainthearted. We'll see how it goes. Interestingly, the order
//  of the base classes (QObject then ncm_vessel_properties) makes a difference
//  in this instance. Also note the keyword 'public' before both parent classes
class ncm_qvessel_properties : public QObject, public ncm_vessel_properties 
{
  Q_OBJECT // needed if we want to use signals or slots

//  INTERFACE

public:
  //  No member variables here, please
  ncm_qvessel_properties();

public slots:
  //: The boolean parameter here is used so that checkboxes and radiobuttons can
  //  interface to these methods via their checked() and toggled() methods
  void setDistal(bool);
  void setNondistal(bool);

  void setSizeUndefined(bool);
  void setSizeNormal(bool);
  void setSizeEnlarged(bool);
  void setSizeGiant(bool);
  void setSizeIrregular(bool);

  void setShapeUndefined(bool);
  void setShapeNormal(bool);

  //: The following set the tortuous/ramified/whatever flag to the value
  //  specified by 'confirm'; all other values are left untouched.
  void setShapeTortuous(bool confirm = true);
  void setShapeRamified(bool confirm = true);

  //: The following set the tortuous/ramified/whatever flag to the value
  //  specified by 'confirm'; all other values are set to zero.
  void setShapeTortuousOnly(bool);
  void setShapeRamifiedOnly(bool);


//  IMPLEMENTATION

protected:

private:

};

#endif // __ncm_qvessel_properties_h__