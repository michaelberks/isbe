#ifndef NAILFOLD_qcapture_end_session_H
#define NAILFOLD_qcapture_end_session_H

#include <QDialog>

#include "ui_ncm_qcapture_end_session.h"
#include <nailfold/qtools/ncm_image_session.h>

class ncm_qcapture_end_session : public QDialog
{
  Q_OBJECT

public:
  ncm_qcapture_end_session(ncm_image_session & session, const QString study_name, const QString subject_study_name,
												QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qcapture_end_session();

public slots:

protected:

private slots:

	void on_sessionNotesPlainTextEdit_textChanged();

	void on_endSessionButton_clicked();
	//void on_


private:
  //  Methods

  //: User interface object
  Ui::end_sessionDialog ui;

	//Current subject (loaded from SQL when subject selected)
	ncm_image_session & session_;
};

#endif // NAILFOLD_qcapture_end_session_H
