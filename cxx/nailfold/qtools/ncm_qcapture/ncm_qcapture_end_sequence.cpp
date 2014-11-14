#include "ncm_qcapture_end_sequence.h"

// VXL libraries
#include <vcl_iostream.h>

// Qt libraries
#include <QDir>
#include <QMessageBox>
#include <QColorDialog>
#include <QtSql/QSqlRecord>
//
//: Constructor
ncm_qcapture_end_sequence::ncm_qcapture_end_sequence(ncm_image_sequence & sequence, const QString study_name, const QString subject_study_name,
																					QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags),
	sequence_(sequence)
{

  // setup the UI
  ui.setupUi(this);
	createHandRadioGroup();

	//Connect up hand radio group
  QObject::connect( handRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(set_digit(int)) );

	//Show study name and initial subject name in UI (the latter can be modified)
	ui.studyLabel->setText("Study: " + study_name);
	ui.nameLabel->setText("Name: " + subject_study_name);
	ui.sequenceNameLineEdit->setText(sequence_.images_sub_dir());
	ui.qualityCheckBox->setChecked(sequence_.acceptable_quality());
}

//
//: Destructor
ncm_qcapture_end_sequence::~ncm_qcapture_end_sequence()
{
}

//
//  Public slots
//

//Update subject's notes
void ncm_qcapture_end_sequence::on_sequenceNotesPlainTextEdit_textChanged()
{
	sequence_.set_notes(ui.sequenceNotesPlainTextEdit->document()->toPlainText());
}
void ncm_qcapture_end_sequence::on_qualityCheckBox_toggled(bool checked)
{
	sequence_.set_acceptable_quality(checked);
}

void ncm_qcapture_end_sequence::on_sequenceNameLineEdit_textChanged()
{
	//Do nothing...
}

//Update subjects dominant hand
void ncm_qcapture_end_sequence::set_digit(int index)
{
	ncm_hand::Hand hand; 
	ncm_hand::Digit digit;
	if (index < 0)
	{
		hand = ncm_hand::NotKnown;
		digit = ncm_hand::None;
	}
	else if (index <= 5)
	{
		hand = ncm_hand::Left;
		digit = ncm_hand::Digit(index);
	}
	else
	{
		hand = ncm_hand::Right;
		digit = ncm_hand::Digit(index - 5);
	}
	sequence_.set_hand(hand);
	sequence_.set_digit(digit);

	//Prefix the sequence name
	QString seq_name = QString(ncm_hand::toLetter(hand).c_str()) + QString(ncm_hand::toLetter(digit).c_str()) + "_"
		+ sequence_.images_sub_dir();
	ui.sequenceNameLineEdit->setText(seq_name);
}

//Return control back to the subject search UI, now with the subject object complete
void ncm_qcapture_end_sequence::on_endSequenceButton_clicked()
{
	if (sequence_.hand() == ncm_hand::NotKnown)
	{
		QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText("No digit been selected");
    msgBox.setInformativeText("Is this OK?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);

    const int response = msgBox.exec();

    if (response == QMessageBox::No)
    {
			return;
		}
	}

	//Check whether a directory with this sequence name exists
	sequence_.set_sequence_name(ui.sequenceNameLineEdit->text());
	if (sequence_.sequence_name() != sequence_.images_sub_dir())
	{
		QDir dir(sequence_.images_root_dir() + "/" + sequence_.sequence_name());
		if (dir.exists())	
		{
			QMessageBox msgBox;
			msgBox.setIcon(QMessageBox::Warning);
			msgBox.setText("Sequence name already exists in this session \n"
				"Please choose a new name");
			msgBox.exec();
			return;
		}
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
void ncm_qcapture_end_sequence::createHandRadioGroup()
{
	handRadioGroup_ = new QButtonGroup(this);
	handRadioGroup_->addButton(ui.otherDigitRadioButton, -2);
	handRadioGroup_->addButton(ui.LD1radioButton, ncm_hand::D1);
	handRadioGroup_->addButton(ui.LD2radioButton, ncm_hand::D2);
	handRadioGroup_->addButton(ui.LD3radioButton, ncm_hand::D3);
	handRadioGroup_->addButton(ui.LD4radioButton, ncm_hand::D4);
	handRadioGroup_->addButton(ui.LD5radioButton, ncm_hand::D5);
	handRadioGroup_->addButton(ui.RD1radioButton, ncm_hand::D1 + 5);
	handRadioGroup_->addButton(ui.RD2radioButton, ncm_hand::D2 + 5);
	handRadioGroup_->addButton(ui.RD3radioButton, ncm_hand::D3 + 5);
	handRadioGroup_->addButton(ui.RD4radioButton, ncm_hand::D4 + 5);
	handRadioGroup_->addButton(ui.RD5radioButton, ncm_hand::D5 + 5);


	//Make exclusive and set 'Not known' as the default checked option
	handRadioGroup_->setExclusive(true);
  ui.otherDigitRadioButton->setChecked(true);

	
}