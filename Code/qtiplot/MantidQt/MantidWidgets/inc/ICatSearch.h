
#ifndef MANTIDWIDGETS_ICATSEARCH_H_
#define MANTIDWIDGETS_ICATSEARCH_H_

//----------------------
// Includes
//----------------------

#include "MantidQtMantidWidgets/ui_ICatSearch.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidQtMantidWidgets/ICatInvestigation.h"
#include "MantidQtMantidWidgets/ICatUtils.h"
#include "WidgetDllOption.h"

#include <algorithm>
#include <QObject>
#include <QCalendarWidget>
#include <QHash>

namespace MantidQt
{
namespace MantidWidgets
{
class  EXPORT_OPT_MANTIDQT_MANTIDWIDGETS ICatSearch : public QWidget
{
  Q_OBJECT

public:
  /// Default Constructor
  ICatSearch(QWidget *parent = 0);
    ~ICatSearch();
 signals:
  void error(const QString&,int param=0);
private:
  /// Initialize the layout
  virtual void initLayout();
  ///populating the instrumentbox 
  void populateInstrumentBox();
  ///get start and end run numbers
  void getRunValues(double& startRun,double& endRun);
  ///get the user set start and end dates
  void getDates(QString& startDate,QString& endDate);
   /// get the user selected instrument
  void getSelectedInstrument(QString& instrName);
  /// execute the algorithm
  bool executeSearch(Mantid::API::ITableWorkspace_sptr& ws1_sptr);

  /// this method populates the search result widget.
  void updatesearchResults(Mantid::API::ITableWorkspace_sptr & ws_sptr );
  /// save settings to registry
  void saveSettings();
  /// read settings from registry
  void readSettings();

  void closeEvent(QCloseEvent* event);

  void setparentWidget(QWidget*);
  //bool eventFilter(QObject *obj, QEvent *event);
  QWidget* getParentWidget();

  //if  casesensitive check box selected 
  bool isCaseSensitiveSearch();

  //event filtering
  bool eventFilter(QObject *obj, QEvent *event);

   /// This method sets the property of algorithm
   template<typename T >
   bool setProperty(QString name,T value)
   {
	
	   try
	   {
		   m_alg->setProperty(name.toStdString(),value);
	   }
	   catch(std::invalid_argument& e)
	   {	
		 
		   emit error(e.what());
           showInvalidMarkerLabel(name);
		 
		   return false;
	   }
	   catch (Mantid::Kernel::Exception::NotFoundError& e)
	   {
		   emit error(e.what());
		   showInvalidMarkerLabel(name);
		   return false;
	   }
	  
	   hideInvalidMarkerLabel(name);
	   return true;
   }

   /// This method adds invalid marker labels to hashtable
   void addtoPropertyLabelsHash();
   /// this method creates shared pointer to search algorithm
   Mantid::API::IAlgorithm_sptr createAlgorithm();
   /// show invalid marker labels
   void showInvalidMarkerLabel(const QString& name);
   /// hide invalid marker labels
   void hideInvalidMarkerLabel(const QString& name);
 
private slots:
	///handler for search button
	void onSearch();
	///handler for close button
	void onClose();
	///
	void investigationSelected(QTableWidgetItem *);
	///start date changed
	void getDate(const QDate& date  );
	///popup DateTime calender to select date
	void popupCalendar();
	//handler for helpbutton
	void helpButtonClicked();
	 
private:
  ///The form generated by Qt Designer
  Ui::ICatSearch m_uiForm;
 
  ///investigation widget
  MantidQt::MantidWidgets::ICatInvestigation* m_invstWidget;

  ///parent widget
  QWidget* m_applicationWindow;

  ///pointer to object to identify starta nd end date tool button
  QObject* m_sender;

    ///stores investigation data 
  Mantid::API::ITableWorkspace_sptr m_ws_sptr;

  ///shared pointer to icat utils object
  boost::shared_ptr<ICatUtils> m_utils_sptr;

  /// hash containing property name and invalid marker * label
  QHash<QString,QWidget*> m_propLabelHash;

  /// shared pointer to search algorithm
  Mantid::API::IAlgorithm_sptr m_alg;


};

}
}

#endif //MANTIDQTCUSTOMINTERFACES_ICATSEARCH_H_
