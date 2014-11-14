#ifndef NAILFOLD_QFINAL_SEARCH_H
#define NAILFOLD_QFINAL_SEARCH_H

#include <QDialog>
#include <QDoubleValidator>

#include "ui_ncm_qfinal_search.h"

class ncm_qfinal_search : public QDialog
{
  Q_OBJECT

public:
  ncm_qfinal_search(QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qfinal_search();

public slots:

protected:

private slots:

private:
  //  Methods

  //: User interface object
  Ui::searchDialog ui;
};

#endif // NAILFOLD_QFINAL_SEARCH_H
