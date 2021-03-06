#include "Utils.h"

namespace swe {

    void printError(const string &message) {

        cout << endl;
        cerr << endl << "[Error] " << message << endl << endl;
        exit(EXIT_FAILURE);
    }

    void printDefault(const string &optionToFind, const string &defaultOption) {

        cout << "[Caution] \n";
        cout << optionToFind + " not specified. Selected DEFAULT: " << defaultOption << ". Please check carefully"
             << " the configuration file if you want to set another " << optionToFind << endl << endl;
    }

    void setupFolder(const string& path, folderOptions_t option) {

        struct stat folderStatus{};
        struct dirent *file;
        vector<string> paths;
        string filepath;
        DIR *folderPtr = opendir(path.c_str());

        // --- Create folder if does not exist, otherwise remove or keep files --- //
        // --- First condition returns 0 whether information about path.c_str() are successfully found --- //
        // --- Second condition returns true whether the directory exists --- //
        if (!(stat(path.c_str(), &folderStatus) == 0 && S_ISDIR(folderStatus.st_mode))) {
            cout << " Creating folder " << path << endl;
            paths = explodePath(path);

            for (auto &i : paths) {
                // --- jump "." and ".." in paths --- //
                if (i != "./" && i != "../")
                    mkdir(i.c_str(), 0755);
            }
        } else {
            if (option == CLEAN_FOLDER) {
                cout << " Cleaning folder " << path << endl;

                while ((file = readdir(folderPtr))) {
                    // --- Build path for each file in the folder --- //
                    filepath = path + "/" + file->d_name;
                    remove(filepath.c_str());
                }

                closedir(folderPtr);
            }
        }
    }

    bool existsInVectorString(const vector<string> &vect, const string &strTarget) {

        if (find(vect.begin(), vect.end(), strTarget) != vect.end())
            return true;
        else
            return false;
    }

    string jumpTo(ifstream &file, const string &keyword, bool &isFound) {

        string line, option;
        isFound = false;

        while (!isFound && !file.eof()) {
            getline(file, line);
            option = line.substr(0, keyword.size());        // --- Extract keyword from line --- //
            if (option == keyword)
                isFound = true;
        }

        return line;
    }

    string trim(const string &str) {

        size_t first = str.find_first_not_of(' ');

        if (string::npos == first)
            return str;

        size_t last = str.find_last_not_of(' ');

        return str.substr(first, (last - first + 1));
    }

    string getOptionString(ifstream &file, const string &optionToFind, vector<string> list, string defaultOption) {

        string line, value;
        bool success;

        // --- Search optionToFind from the beginning of the file --- //
        file.clear();
        file.seekg(0, ios::beg);
        line = jumpTo(file, optionToFind, success);

        // --- Get the selected option --- //
        if (success) {

            value = line.substr(optionToFind.size() + 1, line.size());
            value = trim(value);

            // --- Wrong option typing --- //
            if (find(list.begin(), list.end(), value) == list.end()) {
                cout << endl;
                cerr << endl << "[Error] \n"
                     << value + " --> Wrong " + optionToFind + " option. Select among:" << endl;
                for (auto &j : list) {
                    cout << "  " << j << endl;
                }
                exit(EXIT_FAILURE);
            }
        } else {
            value = move(defaultOption);
            printDefault(optionToFind, value);
        }

        return value;
    }

    vector<string> explodePath(const string &path) {

        int i = 0;
        string currentPath;
        vector<string> paths;
        istringstream ss(path);

        while (getline(ss, currentPath, '/')) {
            currentPath += "/";
            if (i > 0)
                currentPath = paths[i - 1].append(currentPath);
            paths.push_back(currentPath);
            i++;
        }

        return paths;
    }

    vector<string> getOptionVectorStrings(ifstream &file, const string &optionToFind) {

        string line, tmp;
        vector<string> out;
        stringstream ss;
        bool success;

        // --- Search optionToFind from the beginning of the file --- //
        file.clear();
        file.seekg(0, ios::beg);
        line = jumpTo(file, optionToFind, success);

        // --- Get the selected option --- //
        if (success) {
            line = line.substr(optionToFind.size() + 1, line.size());     // --- Read string --- //
            line = trim(line);                                                // --- Trim white spaces --- //
            line = line.substr(1, line.size() - 2);                    // --- Remove brackets --- //
            line = trim(line);                                               // --- Trim white spaces --- //
            ss << line;
            while (getline(ss, tmp, ',')) {
                line = trim(tmp);
                out.push_back(line);
            }
        } else {
            printError("Option " + optionToFind + " not found.");
        }

        return out;
    }

    int getOptionInteger(ifstream &file, const string &optionToFind, bool mandatory, int defaultValue) {

        int value = 0;
        string line, str;
        bool success;

        // --- Search optionToFind from the beginning of the file --- //
        file.clear();
        file.seekg(0, ios::beg);
        line = jumpTo(file, optionToFind, success);

        // --- Get the selected option --- //
        if (success) {
            istringstream ss(line);
            ss >> str >> value;
        } else {
            if (mandatory) {
                printError("Option " + optionToFind + " not found.");
            } else {
                value = defaultValue;
                printDefault(optionToFind, to_string(value));
            }
        }

        return value;
    }

    double getOptionDouble(ifstream &file, const string &optionToFind, bool mandatory, double defaultValue) {

        double value = 0;
        string line, str;
        bool success;

        // --- Search optionToFind from the beginning of the file --- //
        file.clear();
        file.seekg(0, ios::beg);
        line = jumpTo(file, optionToFind, success);

        // --- Get the selected option --- //
        if (success) {
            istringstream ss(line);
            ss >> str >> value;
        } else {
            if (mandatory) {
                printError("Option " + optionToFind + " not found.");
            } else {
                value = defaultValue;
                printDefault(optionToFind, to_string(value));
            }
        }

        return value;
    }
}