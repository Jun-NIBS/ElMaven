#ifndef PROJECTSAVEWORKER_H
#define PROJECTSAVEWORKER_H

#include <QThread>

class MainWindow;
class PeakGroup;

class ProjectSaveWorker : public QThread
{
    Q_OBJECT

public:
    ProjectSaveWorker(MainWindow* mw);

    /**
     * @brief Issue save command for user-specified emDB project.
     * @param fileName Path of the file where the emDB file will be saved.
     * @param saveRawData Whether the file should contain raw EIC and spectra
     * information for peaks.
     */
    virtual void saveProject(const QString fileName,
                             const bool saveRawData = false);

    /**
     * @brief Update the currently set emDB project. This should be the last
     * project that was saved using this thread.
     * @param groupsToSave If any new groups need to be added or modified in
     * this file, they can be sent in this list.
     */
    void updateProject(QList<PeakGroup*> groupsToSave = {});

    /**
     * @brief Obtain the project name for the file that is set as the current
     * project of this thread. Any save commands will save to this file.
     * @return
     */
    QString currentProjectName() const;

protected:
    MainWindow* _mw;
    QString _currentProjectName;
    QList<PeakGroup*> _groupsToSave;
    bool _saveRawData;

    void run();

private:
    /**
     * @brief Write project data into given file as a SQLite database.
     * @details If the file is already a project DB, all project tables are
     * deleted and created anew.
     * @param filename Name of the file to be saved as SQLite project on disk.
     */
    void _saveSqliteProject(const QString fileName);

    /**
     * @brief Save or update the information of a peak group in the current
     * emDB project.
     * @param group Pointer to the `PeakGroup` object which will be saved, or
     * updated.
     * @param filename A string path for filename of the SQLite project to be
     * created (only if it does not exist already).
     */
    void _savePeakGroupInSqlite(PeakGroup *group, QString fileName);
};

class TempProjectSaveWorker : public ProjectSaveWorker
{
    Q_OBJECT

public:
    TempProjectSaveWorker(MainWindow* mw);
    ~TempProjectSaveWorker();

    /**
     * @brief Overrides `ProjectSaveWorker` class' original `saveProject`
     * method, such that calling this will just call the overloaded alternative,
     * thus completely ignoring the file path provided.
     */
    void saveProject(const QString filename,
                     const bool saveRawData = false) override
    {
        saveProject(saveRawData);
    }

    /**
     * @brief This method will save a temporary (time-stamped) emDB file for the
     * current session. If a temporary file has already being saved to, then it
     * will continue to be used until `deleteCurrentProject` is called.
     * @param saveRawData Whether the file should contain raw EIC and spectra
     * information for peaks.
     */
    void saveProject(const bool saveRawData = false);

    /**
     * @brief Deletes the current project (most like the project last saved
     * using this thread) from the filesystem. The next time `saveProject` is
     * called, a new temporary file will be created and save to.
     */
    void deleteCurrentProject();
};

#endif // PROJECTSAVEWORKER_H
