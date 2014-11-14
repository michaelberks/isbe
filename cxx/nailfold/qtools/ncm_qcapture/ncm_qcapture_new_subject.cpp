#include "ncm_qcapture_new_subject.h"

// VXL libraries
#include <vcl_iostream.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_vessel.h>
#include <nailfold/qtools/ncm_qvesselitem.h>

// Qt libraries
#include <QMessageBox>
#include <QColorDialog>
#include <QtSql/QSqlRecord>
//
//: Constructor
ncm_qcapture_new_subject::ncm_qcapture_new_subject(ncm_subject & subject, const int study_id, const QString study_name, QString subject_name,
																					QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags),
	study_id_(study_id),
	study_name_(study_name),
	subject_(subject),
  subject_name_(subject_name)
{
  // setup the UI
  ui.setupUi(this);
	createHandRadioGroup();

	//Connect up hand radio group
  QObject::connect( handRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(set_dominant_hand(int)) );

	//Initialize subject
	subject_.set_study_id(study_id_);
	subject_.set_name_in_study(subject_name_);

	//Show study name and initial subject name in UI (the latter can be modified)
	ui.studyLabel->setText("Study: " + study_name_);
	ui.subjectNameEdit->setText(subject_name_);

}

//
//: Destructor
ncm_qcapture_new_subject::~ncm_qcapture_new_subject()
{
}

//
//  Public slots
//

//
//: Update subject name when text edit changed
void ncm_qcapture_new_subject::on_subjectNameEdit_textChanged()
{
	subject_name_ = ui.subjectNameEdit->text();
	subject_.set_name_in_study(subject_name_);
}

//Update subject's notes
void ncm_qcapture_new_subject::on_subjectNotesPlainTextEdit_textChanged()
{
	subject_.set_notes(ui.subjectNotesPlainTextEdit->document()->toPlainText());
}

//Update subjects dominant hand
void ncm_qcapture_new_subject::set_dominant_hand(int hand)
{
	subject_.set_dominant_hand(ncm_hand::Hand(hand));
}

//Return control back to the subject search UI, now with the subject object complete
void ncm_qcapture_new_subject::on_createSubjectButton_clicked()
{
	if (subject_.dominant_hand() == ncm_hand::NotKnown)
	{
		 QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText("No dominant has been selected");
    msgBox.setInformativeText("Is this OK?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);

    const int response = msgBox.exec();

    if (response == QMessageBox::No)
    {
			return;
		}
	}
	if (subject_.name_in_study().toStdString().empty())
	{
		// Issue a warning if connect failed
      QMessageBox msgBox;
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.setText("Subject name must not be empty.");
      msgBox.exec();
			return;
	}

	//Otherwise we're happy - lose dialog with return code 0 (success)
	done(0);

}
//
//  Private slots
//

//
// Private functions
//
void ncm_qcapture_new_subject::createHandRadioGroup()
{
	handRadioGroup_ = new QButtonGroup(this);
	handRadioGroup_->addButton(ui.nhandRadioButton, ncm_hand::NotKnown);
	handRadioGroup_->addButton(ui.lhandRadioButton, ncm_hand::Left);
	handRadioGroup_->addButton(ui.rhandRadioButton, ncm_hand::Right);

	//Make exclusive and set 'Not known' as the default checked option
	handRadioGroup_->setExclusive(true);
  ui.nhandRadioButton->setChecked(true);

	
}