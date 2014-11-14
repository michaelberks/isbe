#ifndef NAILFOLD_qcapture_new_subject_H
#define NAILFOLD_qcapture_new_subject_H

#include <QDialog>

#include "ui_ncm_qcapture_new_subject.h"
#include <nailfold/qtools/ncm_subject.h>

class ncm_qcapture_new_subject : public QDialog
{
  Q_OBJECT

public:
  ncm_qcapture_new_subject(ncm_subject & subject, const int study_id, const QString study_name, QString subject_study_name,
												QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qcapture_new_subject();

public slots:

protected:

private slots:
	void on_subjectNameEdit_textChanged();
	void set_dominant_hand(int hand);
	void on_subjectNotesPlainTextEdit_textChanged();

	void on_createSubjectButton_clicked();
	//void on_


private:
  //  Methods
	void createHandRadioGroup();

  //: User interface object
  Ui::new_subjectDialog ui;

	//: Radio group to select dominant hand
	QButtonGroup *handRadioGroup_;

	//Currently selected study
	const QString study_name_;
	const int study_id_;

	//Current subject name (in text edit box)
	QString subject_name_;

	//Current subject (loaded from SQL when subject selected)
	ncm_subject & subject_;
};

#endif // NAILFOLD_qcapture_new_subject_H
