#include "ncm_qfinal_search.h"

// VXL libraries
#include <vcl_iostream.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_vessel.h>
#include <nailfold/qtools/ncm_qvesselitem.h>

// Qt libraries
#include <QColorDialog>

//
//: Constructor
ncm_qfinal_search::ncm_qfinal_search(QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags)
{
  // setup the UI
  ui.setupUi(this);
}

//
//: Destructor
ncm_qfinal_search::~ncm_qfinal_search()
{
}

//
//  Public slots
//

//
//: What to do when the user hits 'OK'

//
//  Private slots
//


//
//  Private functions
//
