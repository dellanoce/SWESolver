#ifndef SWE_UTILS_H
#define SWE_UTILS_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>
#include <mpi.h>

#include "Enums.h"

using namespace std;

/*!
 * \brief Enumerator listing options for creating a new directory
 */
enum folderOptions_t {
    PRESERVE_FILES,
    CLEAN_FOLDER
};

/*!
 * \namespace swe
 * \brief Namespace swe contains a list of methods used for parsing the configuration file, for handling errors
 * and for creating new directories.
 * The implementation is based on what is suggested on the following sites:
 * https://stackoverflow.com/questions/8075475/listing-the-files-in-a-directory-and-delete-them-in-c-c
 * https://stackoverflow.com/questions/4980815/c-determining-if-directory-not-a-file-exists-in-linux
 * https://linux.die.net/man/2/stat
 * https://stackoverflow.com/questions/25829143/trim-whitespace-from-a-string/25829178
 * \author Alberto Della Noce
 */

namespace swe {

    /*!
     * \brief Prints an error message and exit the program
     * @param[in] errorMessage - Message to be printed
     */
    void printError(const string &errorMessage);

    /*!
     * \brief Prints a message when a default option is used
     * @param[in] optionToFind - Option not specified
     * @param[in] defaultOption - Default option used
     */
    void printDefault(const string &optionToFind, const string &defaultOption);

    /*!
     * \brief Creates an empty directory if it does not exist. Otherwise it cleans or keeps its files
     * @param[in] path - Folder path
     * @param[in] folderOptions_t - Option for creating directory
     */
    void setupFolder(const string& path, folderOptions_t);

    /*!
     * \brief Returns true whether string strTarget is found in vector vect
     * @param[in] vect - Vector of strings
     * @param[in] strTarget - String to be found
     * @return True or false, depending whether strTarget is found or not
     */
    bool existsInVectorString(const vector<string> &vect, const string &strTarget);

    /*!
     * \brief Returns the line of fileName in which keyword is found
     * @param[in] fileName - Name of the file
     * @param[in] keyword - Word to be searched
     * @param[out] isFound - True or false, depending whether keyword is found or not
     * @return Line containing keyword
     */
    string jumpTo(ifstream &fileName, const string &keyword, bool &isFound);

    /*!
     * \brief Trims spaces in a string
     * @param[in] str - Input string
     * @return String with no spaces
     */
    string trim(const string &str);

    /*!
     * \brief Returns the selected string option among the possible contained in the list
     * @param[in] fileName - Name of the file
     * @param[in] optionToFind - Option to be found
     * @param[in] list - List of possible options
     * @param[in] defaultOption - Option set by default
     * @return Selected option. Otherwise returns the default option
     */
    string getOptionString(ifstream &fileName, const string &optionToFind, vector<string> list, string defaultOption);

    /*!
     * \brief Returns a vector containing all path parts
     * @param[in] path - Path to be exploded
     * @return Vector containing path sub-paths
     */
    vector<string> explodePath(const string &path);

    /*!
     * \brief Returns the vector of strings option
     * @param[in] fileName - Name of the file
     * @param[in] optionToFind - Option to be found
     * @return Vector containing values of optionToFind
     */
    vector<string> getOptionVectorStrings(ifstream &fileName, const string &optionToFind);

    /*!
     * \brief Returns the selected integer option. If it is mandatory, returns the defaultValue if the option
     * is not found
     * @param[in] fileName - Name of the file
     * @param[in] optionToFind - Option to be found
     * @param[in] mandatory - Bool stating whether the option is mandatory or not
     * @param[in] defaultValue - Default value
     * @return Selected option. Otherwise returns the default option
     */
    int getOptionInteger(ifstream &fileName, const string &optionToFind, bool mandatory, int defaultValue = 0);

    /*!
     * \brief Returns the selected double option. If it is mandatory, returns the defaultValue if the option
     * is not found
     * @param[in] fileName - Name of the file
     * @param[in] optionToFind - Option to be found
     * @param[in] mandatory - Bool stating whether the option is mandatory or not
     * @param[in] defaultValue - Default value
     * @return Selected option. Otherwise returns the default option
     */
    double getOptionDouble(ifstream &fileName, const string &optionToFind, bool mandatory, double defaultValue = 0.0);

// ----------------------------------------------- TEMPLATES ------------------------------------------------- //

    /*!
     * \brief Returns the enumerator related to the option string within the given map
     * @tparam[in] T - Datatype for sub-options enumerator
     * @param[in] option - Option to be found
     * @param[in] optionsMap - Map containing sub-options
     * @return Sub-option enumerator
     */
    template<typename T>
    T getEnumFromMap(string option, map<T, string> optionsMap) {

        bool isFound = false;
        T enumToSet;

        for (auto const &elem : optionsMap) {
            if (elem.second == option) {
                enumToSet = elem.first;
                isFound = true;
            }
        }

        if (!isFound)
            printError("Option " + option + " does not match any enum identifier.");

        return enumToSet;
    }
}

#endif //SWE_UTILS_H
