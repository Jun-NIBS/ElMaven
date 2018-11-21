#ifndef PROJECTDOCKWIDGET_H
#define PROJECTDOCKWIDGET_H

#include "stable.h"
#include "mzSample.h"
#include "mainwindow.h"
#include <boost/signals2.hpp>
#include <boost/bind.hpp>

extern Database DB; 

class ProjectDockWidget : public QDockWidget
{
    Q_OBJECT
public:
    explicit ProjectDockWidget(QMainWindow *parent = 0);
    boost::signals2::signal< void (const string&,unsigned int , int ) > boostSignal;
    QString lastOpenedProject;
    QString lastSavedProject;
    QColor  lastUsedSampleColor;
    QMap<mzSample*, QColor> storeSampleColors;
    QTreeWidget* getTreeWidget();
    void prepareSampleCohortFile(QString sampleCohortFileName);

Q_SIGNALS:

public Q_SLOTS:
    void showInfo();
    QString getProjectDescription();
    void setProjectDescription(QString text);
    void setInfo(vector<mzSample*>&samples);
    void changeSampleOrder();
    void updateSampleList();

    /**
     * @brief Bring up a file dialog and allow the user to save current session
     * as an emDB project.
     */
    void saveProjectAsSQLite();

    /**
     * @brief Update saved data in the currently open emDB project.
     */
    void saveSQLiteProject();

    /**
     * @brief Save a bookmarked PeakGroup in currently open emDB project.
     * @param group The bookmarked PeakGroup to be saved.
     * @return Unique integer ID of the saved bookmark if successful, -1
     * otherwise.
     */
    int saveBookmarkedGroup(PeakGroup* group);

    /**
     * @brief Load session data saved within an emDB SQLite project.
     * @param filename String name of the emDB file to be loaded.
     */
    void loadSQLiteProject(QString filename);

    /**
     * @brief Save any pending changes and close the currently open SQLite
     * (emDB) project.
     */
    void saveAndCloseCurrentSQLiteProject();

    /**
     * @brief Clear out structures before loading a new SQLite (emDB) project.
     * @details This method tries to unload any loaded samples and deletes
     * peak tables in the current session. Although it should be noted that
     * this has historically not resulted in truly fresh session experience for
     * users. Still something like this should ideally be available to
     * developers for house cleaning before loading new projects into existing
     * sessions.
     */
    void clearSession();

    void saveProjectAsMzRoll();
    void saveMzRollProject();
    void loadMzRollProject(QString filename);
    void saveMzRollProject(QString filename, TableDockWidget* peakTable=nullptr);
    void setSampleColor(mzSample* sample, QColor color); //TODO: Sahil, Added while merging projectdockwidget
    void unloadSelectedSamples(); //TODO: Sahil, Added while merging projectdockwidget
    void sendBoostSignal( const string& progressText, unsigned int completed_samples, int total_samples)
    {
    boostSignal(progressText, completed_samples, total_samples);
    }


    void SetAsBlankSamples(); //TODO: --@Giridhari, Create function to Set samples as Blank Samples
protected Q_SLOTS:
      void keyPressEvent( QKeyEvent *e );

private Q_SLOTS:
    void showSample(QTreeWidgetItem* item, int col);
    void showSampleInfo(QTreeWidgetItem* item, int col);
    void changeSampleColor(QTreeWidgetItem* item, int col);
    void changeNormalizationConstant(QTreeWidgetItem* item, int col);
    void changeSampleSet(QTreeWidgetItem* item, int col);
    void selectSample(QTreeWidgetItem* item, int col);
    void changeColors();
    void checkUncheck(); //TODO: Sahil, Added while merging projectdockwidget
    void setSampleColor(QTreeWidgetItem* item, QColor color);
    void dropEvent (QDropEvent*event);
    // void unloadSample(); //TODO: Sahil, Removed while merging projectdockwidget
    void unloadSample(mzSample* sample); //TODO: Sahil, Added while merging projectdockwidget
    void filterTreeItems(QString filterString);

private:
    QTreeWidgetItem* getParentFolder(QString filename);
    void markBlank(QTreeWidgetItem* item);
    void unmarkBlank(QTreeWidgetItem* item);
    QMap<QString,QTreeWidgetItem*> parentMap;
    QTextEdit* _editor;
    MainWindow* _mainwindow;
    QTreeWidget* _treeWidget;


    QMap<QString, QColor> storeColor;
    QColor  usedColor;



};

#endif // PROJECTDOCKWIDGET_H
