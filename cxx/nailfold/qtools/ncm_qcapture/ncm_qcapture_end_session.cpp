#include "ncm_qcapture_end_session.h"

// VXL libraries
#include <vcl_iostream.h>

// Qt libraries
#include <QMessageBox>
#include <QColorDialog>
#include <QtSql/QSqlRecord>
//
//: Constructor
ncm_qcapture_end_session::ncm_qcapture_end_session(ncm_image_session & session, const QString study_name, const QString subject_study_name,
																					QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags),
	session_(session)
{
  // setup the UI
  ui.setupUi(this);

	//Show study name and initial subject name in UI (the latter can be modified)
	ui.studyLabel->setText("Study: " + study_name);
	ui.nameLabel->setText("Name: " + subject_study_name);
	ui.summaryTextEdit->setText(session_.sequences_string());
	ui.summaryTextEdit->setReadOnly(true);
}

//
//: Destructor
ncm_qcapture_end_session::~ncm_qcapture_end_session()
{
}

//
//  Public slots
//

//Update subject's notes
void ncm_qcapture_end_session::on_sessionNotesPlainTextEdit_textChanged()
{
	session_.set_notes(ui.sessionNotesPlainTextEdit->document()->toPlainText());
}

//Return control back to the subject search UI, now with the subject object complete
void ncm_qcapture_end_session::on_endSessionButton_clicked()
{
	//Otherwise we're happy - lose dialog with return code 0 (success)
	done(0);
}
//
//  Private slots
//

//
// Private functions
//