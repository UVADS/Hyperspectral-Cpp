#include "write_log_file.h"
#include "settings.h"
#include <iostream>
#include <fstream>
using namespace std;

string getCurrentDateTime( string s ){
    time_t now = time(0);
    struct tm  tstruct;
    char  buf[80];
    tstruct = *localtime(&now);
    if(s=="now")
        strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
    else if(s=="date")
        strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);
    return string(buf);
}

void Logger(string logMsg, settingsValues settings){
    string now = getCurrentDateTime("now");
    fstream file;
    file.open(settings.loggerPath, std::ios_base::out | std::ios_base::app );
    file << now << '\t' << logMsg << endl;
    file.close();
    if (settings.verbose > 0){
        cout << logMsg << endl;
    }
}

void Logger_initialize(settingsValues settings){
    cout << "LOG FILE NAME: " << settings.loggerPath << endl;
    ofstream ofs;
    ofs.open(settings.loggerPath, std::ofstream::out | std::ofstream::trunc);
    ofs.close();
}