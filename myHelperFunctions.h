// This file contains heplfull functions like reading from a textfile

using namespace std;


//---------------------FORWARD DECLARATIONS----------------------------------------------------------------

void PrintText(string text);
tuple<vector<double>, vector<double> , bool> readMatrixMSR(const std::string& filename);

//------------------------------------------------------------------------------------------------------------------------

void PrintText(string text)
{
    cout << text << endl;
}

//------------------------------------------------------------------------------------------------------------------------

tuple<vector<double>, vector<double>, bool> readMatrixMSR(const std::string& filename) {
    std::ifstream file(filename); // Open the file
    std::vector<double> array1;  //JM
    std::vector<double> array2;  //VM
    bool sym = 0;

    if (file.is_open()) {
        std::string line;
        int count = 0;  //used to acces first line
        while (std::getline(file, line)) {
            if (!line.empty()) {
                std::istringstream iss(line);
                double num1;
                double num2;


                if( count == 0)
                {
                    if (line == " s")
                    {
                        cout << "First line : " << line << " . Hence Matrix is symmetric!" << endl;
                        sym = 1;
                    }
                    count ++;
                }

                if (iss >> num1 >> num2) {
                    if(count == 1)
                    {
                        count++;
                    }
                    else
                    {
                        array1.push_back(num1);
                        array2.push_back(num2);
                    }
                }

            }
        }
        file.close(); // Close the file
    }
    return {array1,array2,sym};
}

//------------------------------------------------------------------------------------------------------------------------
