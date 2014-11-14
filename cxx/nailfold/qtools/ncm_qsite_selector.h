#ifndef ncm_qsite_selector_h
#define ncm_qsite_selector_h

#include <QDialog>
#include "ui_ncm_qsite_selector.h"

class ncm_qsite_selector : public QDialog
{
  Q_OBJECT

public:
  ncm_qsite_selector(QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qsite_selector();

public slots:
  virtual void accept();

protected:

private slots:

private:
  //  Methods
  QButtonGroup* radioGroup_;

  //: User interface object
  Ui::Dialog ui;
};

#endif // ncm_qsite_selector_h
