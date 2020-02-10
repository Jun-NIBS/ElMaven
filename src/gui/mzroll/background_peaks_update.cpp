#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <QDir>
#include <QCoreApplication>
#include <QProcess>
#include <QJsonObject>

#include "Compound.h"
#include "adductdetection.h"
#include "alignmentdialog.h"
#include "common/analytics.h"
#include "csvreports.h"
#include "background_peaks_update.h"
#include "database.h"
#include "grouprtwidget.h"
#include "isotopeDetection.h"
#include "mainwindow.h"
#include "masscutofftype.h"
#include "mavenparameters.h"
#include "mzAligner.h"
#include "mzSample.h"
#include "obiwarp.h"
#include "PeakDetector.h"
#include "samplertwidget.h"

BackgroundPeakUpdate::BackgroundPeakUpdate(QWidget*) {
        mainwindow = NULL;
        _stopped = true;
        setTerminationEnabled(true);
        runFunction = "computeKnowsPeaks";
        peakDetector = nullptr;
        setPeakDetector(new PeakDetector());
        _classifyPeakGroup = true;
}

BackgroundPeakUpdate::~BackgroundPeakUpdate() {
    delete peakDetector;
    mavenParameters->cleanup();  // remove allgroups
}

/**
 * BackgroundPeakUpdate::run This function starts the thread. This function is
 * called by start() internally in QTThread. start() function will be called
 * where the thread starts.
 */
void BackgroundPeakUpdate::run(void) {
        //Making sure that instance of the mainwindow is present so that the
        //peakdetection process can be ran
        if (mainwindow == NULL) {
                quit();
                return;
        }
        connect(this, SIGNAL(alignmentError(QString)), mainwindow, SLOT(showAlignmentErrorDialog(QString)));
        if (mavenParameters->alignSamplesFlag) {
                connect(this, SIGNAL(alignmentComplete(QList<PeakGroup> )), mainwindow, SLOT(showAlignmentWidget()));
        }
	qRegisterMetaType<QList<PeakGroup> >("QList<PeakGroup>");
	connect(this, SIGNAL(alignmentComplete(QList<PeakGroup> )), mainwindow, SLOT(plotAlignmentVizAllGroupGraph(QList<PeakGroup>)));
	connect(this, SIGNAL(alignmentComplete(QList<PeakGroup> )), mainwindow->groupRtWidget, SLOT(setCurrentGroups(QList<PeakGroup>)));
        connect(this, SIGNAL(alignmentComplete(QList<PeakGroup> )), mainwindow->sampleRtWidget, SLOT(plotGraph()));
        mavenParameters->stop = false;
        started();
        if (runFunction == "alignUsingDatabase") {
                alignUsingDatabase();
        } else if (runFunction == "processSlices") {
            processSlices();
        } else if (runFunction == "processMassSlices") {
                processMassSlices();
        } else if (runFunction == "pullIsotopes") {
                pullIsotopes(mavenParameters->_group);
        } else if (runFunction == "pullIsotopesIsoWidget") {
                pullIsotopesIsoWidget(mavenParameters->_group);
        } else if (runFunction == "pullIsotopesBarPlot") {
                pullIsotopesBarPlot(mavenParameters->_group);
        } else if (runFunction == "computePeaks") {
                computePeaks();
        } else if(runFunction == "alignWithObiWarp" ){
                alignWithObiWarp();
        } else {
                qDebug() << "Unknown Function " << runFunction.c_str();
        }

        quit();
        return;
}
void BackgroundPeakUpdate::alignWithObiWarp()
{
    ObiParams *obiParams = new ObiParams(mainwindow->alignmentDialog->scoreObi->currentText().toStdString(),
                                         mainwindow->alignmentDialog->local->isChecked(),
                                         mainwindow->alignmentDialog->factorDiag->value(),
                                         mainwindow->alignmentDialog->factorGap->value(),
                                         mainwindow->alignmentDialog->gapInit->value(),
                                         mainwindow->alignmentDialog->gapExtend->value(),
                                         mainwindow->alignmentDialog->initPenalty->value(),
                                         mainwindow->alignmentDialog->responseObiWarp->value(),
                                         mainwindow->alignmentDialog->noStdNormal->isChecked(),
                                         mainwindow->alignmentDialog->binSizeObiWarp->value());

    Q_EMIT(updateProgressBar("Aligning samples…", 0, 100));

    Aligner aligner;
    aligner.setAlignmentProgress.connect(boost::bind(&BackgroundPeakUpdate::qtSlot,
                                                     this, _1, _2, _3));

    _stopped = aligner.alignWithObiWarp(mavenParameters->samples, obiParams, mavenParameters);
    delete obiParams;

    if (_stopped) {
        Q_EMIT(restoreAlignment());
        //restore previous RTs
        for (auto sample : mavenParameters->samples) {
            sample->restorePreviousRetentionTimes();
        }

        //stopped without user intervention
        if (!mavenParameters->stop)
            Q_EMIT(alignmentError(
                   QString("There was an error during alignment. Please try again.")));

        mavenParameters->stop = false;
        return;
    }
        
    mainwindow->sampleRtWidget->plotGraph();
    Q_EMIT(samplesAligned(true));

}
void BackgroundPeakUpdate::writeCSVRep(string setName)
{
    peakDetector->pullAllIsotopes();

    if (_classifyPeakGroup)
        classifyGroups(mavenParameters->allgroups);

    for (int j = 0; j < mavenParameters->allgroups.size(); j++) {
        if (mavenParameters->keepFoundGroups) {
            Q_EMIT(newPeakGroup(&(mavenParameters->allgroups[j])));
            QCoreApplication::processEvents();
        }
    }
    Q_EMIT(updateProgressBar("Done", 1, 1));
}

void BackgroundPeakUpdate::setPeakDetector(PeakDetector *pd)
{
    if (peakDetector != nullptr)
        delete peakDetector;

    peakDetector = pd;
    peakDetector->boostSignal.connect(boost::bind(&BackgroundPeakUpdate::qtSlot,
                                                  this,
                                                  _1,
                                                  _2,
                                                  _3));
}

void BackgroundPeakUpdate::processSlices() {
        processSlices(mavenParameters->_slices, "sliceset");
}

void BackgroundPeakUpdate::processSlice(mzSlice& slice) {
        vector<mzSlice*> slices;
        slices.push_back(&slice);
        processSlices(slices, "sliceset");
}

//TODO: kiran Make a function which tell that its from option
//window and should be called where the settings are been called
void BackgroundPeakUpdate::getProcessSlicesSettings() {
        QSettings* settings = mainwindow->getSettings();

        // To Do: Are these lines required. The same is already being done in PeakDetectionDialog.cpp
        // mavenParameters->baseline_smoothingWindow = settings->value(
        //         "baseline_smoothing").toInt();
        // mavenParameters->baseline_dropTopX =
        //         settings->value("baseline_quantile").toInt();

}

void BackgroundPeakUpdate::align() {

        //These else if statements will take care of all corner cases of undoAlignment
        if (mavenParameters->alignSamplesFlag && mavenParameters->alignButton > 0) {
                ;
        } else if (mavenParameters->alignSamplesFlag && mavenParameters->alignButton ==0){
                mavenParameters->alignButton++;
                mavenParameters->undoAlignmentGroups = mavenParameters->allgroups;
        } else if (mavenParameters->alignSamplesFlag && mavenParameters->alignButton == -1) {
                ;
        } else {
                mavenParameters->alignButton = -1;
                mavenParameters->undoAlignmentGroups = mavenParameters->allgroups;
        }

        if (mavenParameters->alignSamplesFlag && !mavenParameters->stop) {
                Q_EMIT(updateProgressBar("Aligning samples…", 0, 0));
                vector<PeakGroup*> groups(mavenParameters->allgroups.size());
                for (int i = 0; i < mavenParameters->allgroups.size(); i++)
                        groups[i] = &mavenParameters->allgroups[i];
                Aligner aligner;
                int alignAlgo = mainwindow->alignmentDialog->alignAlgo->currentIndex();

                if (alignAlgo == 1) {
                        mainwindow->getAnalytics()->hitEvent("Alignment", "PolyFit");
                        aligner.setMaxIterations(mainwindow->alignmentDialog->maxIterations->value());
                        aligner.setPolymialDegree(mainwindow->alignmentDialog->polynomialDegree->value());
                        aligner.doAlignment(groups);
                        mainwindow->sampleRtWidget->setDegreeMap(aligner.sampleDegree);
                        mainwindow->sampleRtWidget->setCoefficientMap(aligner.sampleCoefficient);
                }

                mainwindow->deltaRt = aligner.getDeltaRt();
                mavenParameters->alignSamplesFlag = false;

        }
        QList<PeakGroup> listGroups;
        for (unsigned int i = 0; i<mavenParameters->allgroups.size(); i++) {
                listGroups.append(mavenParameters->allgroups.at(i));
        }	

        Q_EMIT(alignmentComplete(listGroups));
        Q_EMIT(samplesAligned(true));
}

void BackgroundPeakUpdate::alignUsingDatabase() {

    vector<mzSlice*> slices =
        peakDetector->processCompounds(mavenParameters->compounds, "compounds");
        processSlices(slices, "compounds");


}

void BackgroundPeakUpdate::processSlices(vector<mzSlice*>&slices,
                                         string setName) {
        getProcessSlicesSettings();

        peakDetector->processSlices(slices, setName);

        if (runFunction == "alignUsingDatabase") align();

        if (mavenParameters->showProgressFlag
            && mavenParameters->pullIsotopesFlag) {
                Q_EMIT(updateProgressBar("Calculating Isotopes", 1, 100));
        }

        if (mavenParameters->searchAdducts) {
            mavenParameters->allgroups = AdductDetection::findAdducts(
                mavenParameters->allgroups,
                mavenParameters->getChosenAdductList(),
                peakDetector);

            // we remove adducts for which parent ion's RT is too different
            AdductDetection::filterAdducts(mavenParameters->allgroups,
                                           mavenParameters);
            emit (updateProgressBar("Filtering out false adducts…", 0, 0));
        }

        writeCSVRep(setName);
}

void BackgroundPeakUpdate::qtSlot(const string& progressText, unsigned int progress, int totalSteps)
{
        Q_EMIT(updateProgressBar(QString::fromStdString(progressText), progress, totalSteps));

}

void BackgroundPeakUpdate::processCompounds(vector<Compound*> set,
                                            string setName)
{
    if (set.size() == 0)
            return;

    Q_EMIT(updateProgressBar("Processing Compounds", 0, 0));

    vector<mzSlice*> slices = peakDetector->processCompounds(set, setName);
    processSlices(slices, setName);
    delete_all(slices);
}

void BackgroundPeakUpdate::processMassSlices() {
        Q_EMIT (updateProgressBar("Computing Mass Slices", 0, 0));
        mavenParameters->sig.connect(boost::bind(&BackgroundPeakUpdate::qtSignalSlot, this, _1, _2, _3));
        peakDetector->processMassSlices(mavenParameters->compounds);

        align();

        writeCSVRep("allslices");
}

void BackgroundPeakUpdate::qtSignalSlot(const string& progressText, unsigned int completed_slices, int total_slices)
{
        Q_EMIT(updateProgressBar(QString::fromStdString(progressText), completed_slices, total_slices));

}

void BackgroundPeakUpdate::completeStop() {

    peakDetector->resetProgressBar();
        mavenParameters->stop = true;
        stop();
}

void BackgroundPeakUpdate::computePeaks() {

        processCompounds(mavenParameters->compounds, "groups");
}

/**
 * BackgroundPeakUpdate::setRunFunction Getting the function that has tobe ran
 * as a thread and updating it inside a variable
 * @param functionName [description]
 */
void BackgroundPeakUpdate::setRunFunction(QString functionName) {
        runFunction = functionName.toStdString();
}

void BackgroundPeakUpdate::pullIsotopes(PeakGroup* parentgroup) {

	bool isotopeFlag = mavenParameters->pullIsotopesFlag;

        if (!isotopeFlag) return;

	bool C13Flag = mavenParameters->C13Labeled_BPE;
	bool N15Flag = mavenParameters->N15Labeled_BPE;
	bool S34Flag = mavenParameters->S34Labeled_BPE;
	bool D2Flag = mavenParameters->D2Labeled_BPE;

        IsotopeDetection::IsotopeDetectionType isoType;
        isoType = IsotopeDetection::PeakDetection;

	IsotopeDetection isotopeDetection(
                mavenParameters,
                isoType,
                C13Flag,
                N15Flag,
                S34Flag,
                D2Flag);
	isotopeDetection.pullIsotopes(parentgroup);
}

void BackgroundPeakUpdate::pullIsotopesIsoWidget(PeakGroup* parentgroup) {
	bool C13Flag = mavenParameters->C13Labeled_BPE;
	bool N15Flag = mavenParameters->N15Labeled_BPE;
	bool S34Flag = mavenParameters->S34Labeled_BPE;
	bool D2Flag = mavenParameters->D2Labeled_BPE;

        IsotopeDetection::IsotopeDetectionType isoType;
        isoType = IsotopeDetection::IsoWidget;

	IsotopeDetection isotopeDetection(
                mavenParameters,
                isoType,
                C13Flag,
                N15Flag,
                S34Flag,
                D2Flag);
	isotopeDetection.pullIsotopes(parentgroup);
}

void BackgroundPeakUpdate::pullIsotopesBarPlot(PeakGroup* parentgroup) {

	bool C13Flag = mavenParameters->C13Labeled_Barplot;
	bool N15Flag = mavenParameters->N15Labeled_Barplot;
	bool S34Flag = mavenParameters->S34Labeled_Barplot;
	bool D2Flag = mavenParameters->D2Labeled_Barplot;

        IsotopeDetection::IsotopeDetectionType isoType;
        isoType = IsotopeDetection::BarPlot;

	IsotopeDetection isotopeDetection(
                mavenParameters,
                isoType,
                C13Flag,
                N15Flag,
                S34Flag,
                D2Flag);
        isotopeDetection.pullIsotopes(parentgroup);
}

bool BackgroundPeakUpdate::covertToMzXML(QString filename, QString outfile) {

        QFile test(outfile);
        if (test.exists())
                return true;

        QString command = QString("ReAdW.exe --centroid --mzXML \"%1\" \"%2\"").arg(
                filename).arg(outfile);

        qDebug() << command;

        QProcess *process = new QProcess();
        //connect(process, SIGNAL(finished(int)), this, SLOT(doVideoCreated(int)));
        process->start(command);

        if (!process->waitForStarted()) {
                process->kill();
                return false;
        }

        while (!process->waitForFinished()) {
        };
        QFile testOut(outfile);
        return testOut.exists();
}

void BackgroundPeakUpdate::classifyGroups(vector<PeakGroup>& groups)
{
    Q_EMIT(updateProgressBar("Classifying peaks…", 0, 0));
    if (groups.empty())
        return;

    auto tempDir = QStandardPaths::writableLocation(
                       QStandardPaths::GenericConfigLocation)
                   + QDir::separator()
                   + "ElMaven";

    // TODO: binary name will keep changing and should not be hardcoded
    // TODO: model should not exist anywhere on the filesystem
#if defined(Q_OS_MAC) || defined(Q_OS_LINUX)
    auto mlBinary = tempDir + QDir::separator() + "moi";
#endif
#ifdef Q_OS_WIN
    auto mlBinary = tempDir + QDir::separator() + "moi.exe";
#endif
    auto mlModel = tempDir + QDir::separator() + "model.pickle.dat";

    if (!QFile::exists(mlBinary)) {
        cerr << "Error: ML binary not found at path: "
             << mlBinary.toStdString()
             << endl;
        return;
    }
    if (!QFile::exists(mlModel)) {
        cerr << "Error: ML model not found at path: "
             << mlModel.toStdString()
             << endl;
        return;
    }

    QString peakAttributesFile = tempDir
                                 + QDir::separator()
                                 + "peak_ml_input.csv";
    QString classificationOutputFile = tempDir
                                       + QDir::separator()
                                       + "peak_ml_output.csv";

    // have to enumerated and assign each group with an ID, because multiple
    // groups at this point have the same ID
    int startId = 1;
    for (auto& group : groups)
        group.groupId = startId++;
    CSVReports::writeDataForPeakMl(peakAttributesFile.toStdString(),
                                   groups);
    if (!QFile::exists(peakAttributesFile)) {
        cerr << "Error: peak attributes input file not found at path: "
             << peakAttributesFile.toStdString()
             << endl;
        return;
    }

    QStringList mlArguments;
    mlArguments << "--input_attributes_file" << peakAttributesFile
                << "--output_moi_file" << classificationOutputFile
                << "--model_path" << mlModel;
    QProcess subProcess;
    subProcess.setWorkingDirectory(QFileInfo(mlBinary).path());
    subProcess.start(mlBinary, mlArguments);
    subProcess.waitForFinished(-1);
    if (!QFile::exists(classificationOutputFile)) {
        qDebug() << "MOI stdout:" << subProcess.readAllStandardOutput();
        qDebug() << "MOI stderr:" << subProcess.readAllStandardError();
        cerr << "Error: peak classification output file not found at path: "
             << peakAttributesFile.toStdString()
             << endl;
        return;
    }

    QFile file(classificationOutputFile);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        cerr << "Error: failed to open classification output file" << endl;
        return;
    }

    map<int, pair<int, float>> predictions;
    map<int, multimap<float, string>> inferences;
    map<int, vector<pair<int, float>>> correlations;
    map<string, int> headerColumnMap;
    map<string, int> attributeHeaderColumnMap;
    vector<string> headers;
    vector<string> knownHeaders = {"groupId",
                                   "label",
                                   "probability",
                                   "correlations"};
    while (!file.atEnd()) {
        string line = file.readLine().trimmed().toStdString();
        if (line.empty())
            continue;

        vector<string> fields;
        mzUtils::splitNew(line, ",", fields);
        for (auto& field : fields) {
            int n = field.length();
            if (n > 2 && field[0] == '"' && field[n-1] == '"') {
                field = field.substr(1, n-2);
            }
            if (n > 2 && field[0] == '\'' && field[n-1] == '\'') {
                field = field.substr(1, n-2);
            }
        }

        if (headers.empty()) {
            headers = fields;
            for(size_t i = 0; i < fields.size(); ++i) {
                if (find(begin(knownHeaders),
                         end(knownHeaders),
                         fields[i]) != knownHeaders.end()) {
                    headerColumnMap[fields[i]] = i;
                } else {
                    attributeHeaderColumnMap[fields[i]] = i;
                }
            }
            continue;
        }

        int groupId = -1;
        int label = -1;
        float probability = 0.0f;
        vector<pair<int, float>> group_correlations;
        if (headerColumnMap.count("groupId"))
            groupId = string2integer(fields[headerColumnMap["groupId"]]);
        if (headerColumnMap.count("label"))
            label = string2integer(fields[headerColumnMap["label"]]);
        if (headerColumnMap.count("probability"))
            probability = string2float(fields[headerColumnMap["probability"]]);
        if (headerColumnMap.count("correlations")) {
            QString correlations_str =
                QString::fromStdString(fields[headerColumnMap["correlations"]]);
            if (correlations_str != "[]") {
                // trim brackets from start and end
                auto size = correlations_str.size();
                correlations_str = correlations_str.mid(2, size - 4);

                // obtain individual {group ID <-> correlation score} pair
                QStringList pair_strings = correlations_str.split(") (");
                for (auto& pair_str : pair_strings) {
                    QStringList pair = pair_str.split(" ");
                    group_correlations.push_back(
                        make_pair(string2integer(pair.at(0).toStdString()),
                                  string2float(pair.at(1).toStdString())));
                }
            }
        }

        if (groupId != -1) {
            predictions[groupId] = make_pair(label, probability);

            // we use multimap for values mapping to attribute names because
            // later on sorted values are really helpful
            multimap<float, string> groupInference;
            for (auto& element : attributeHeaderColumnMap) {
                string attribute = element.first;
                float value = string2float(fields[element.second]);
                groupInference.insert(make_pair(value, attribute));
            }
            inferences[groupId] = groupInference;
            correlations[groupId] = group_correlations;
        }
    }

    for (auto& group : groups) {
        if (predictions.count(group.groupId)) {
            pair<int, float> prediction = predictions.at(group.groupId);
            group.setPredictedLabel(
                PeakGroup::classificationLabelForValue(prediction.first),
                prediction.second);
        }
        if (inferences.count(group.groupId))
            group.setPredictionInference(inferences.at(group.groupId));

        // add correlated groups
        auto& group_correlations = correlations.at(group.groupId);
        for (auto& elem : group_correlations) {
            auto groupId = elem.first;
            auto correlationFactor = elem.second;
            auto result = find_if(begin(groups),
                                  end(groups),
                                  [groupId] (const PeakGroup& other) {
                                      return other.groupId == groupId;
                                  });
            if (result != end(groups))
                group.addCorrelatedGroup(&(*result), correlationFactor);
        }
    }
}
