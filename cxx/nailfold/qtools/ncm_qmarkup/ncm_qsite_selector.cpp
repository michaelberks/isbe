#include "ncm_qsite_selector.h"

// VXL libraries
#include <vcl_iostream.h>
#include <vcl_string.h>

//
//: Constructor
ncm_qsite_selector::ncm_qsite_selector(QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags)
{
  // setup the UI
  ui.setupUi(this);

  radioGroup_ = new QButtonGroup;
  for (int i = 0; i < ui.groupBox->children().size(); ++i)
  {
    QRadioButton* radio = 
        qobject_cast<QRadioButton*>(ui.groupBox->children()[i]);

    if (radio != NULL)
      radioGroup_->addButton(radio);
  }
  radioGroup_->setExclusive(true);
}

//
//: Destructor
ncm_qsite_selector::~ncm_qsite_selector()
{
}

//
//: What to do when the user hits 'OK'
void ncm_qsite_selector::accept()
{
  // Radiobuttons have names in format <initial><surname>Radio.
  // Therefore, chop off last five characters to get a valid username.
  QString username = radioGroup_->checkedButton()->objectName();
  username.chop(5);

  // Write site name to a text file
  vcl_ofstream ofs("./site.txt");
  ofs << username.toStdString() << '\n';
  ofs.close();

  // close dialog with return code 0 (success)
  done(0);
}


