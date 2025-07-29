// This file contains all needed mathematicall functions
// like for exaple vector math like dotproduct and vector addition.
// there are also functions to cut vectors
using namespace std;

//--------------------Forwar declarations--------------------------------------------------
vector<double> MatVecProd(bool& sym, vector <double>& JM , vector <double>& VM, vector<double>& x);
vector<double> residualCalc(bool& sym,vector<double>& JM, vector<double>& VM, vector<double>& x, vector<double>& b);
void PrintVec(vector<double>& vec);
vector<double> VecAdd(vector<double> vec1,vector<double> vec2);
vector<double> VecSub(vector<double> vec1, vector<double> vec2);
vector<double> PickSlice(vector<double> vec,vector<double> indices);
double DotProd(vector<double> vec1, vector<double> vec2);
vector<double> VecSlice(vector<double>& vec,int start, int end);
double VecNorm(vector<double> vec);
vector<double> VecElemDiv(vector<double> vec, double value);
vector<double> VecNormalize(vector<double> vec);
void PrintMatrix(vector<vector<double>> mat);
vector<vector<double>> CreateMatrix(double row, double coloumn);
vector<double> VecScalarMult(vector<double> vec , double scalar);
vector<vector<double>> MatrixTrans(vector<vector<double>>& mat);
vector<double> backSub(vector<vector<double>>& mat , vector<double>& b);
vector<vector<double>> MatrixTruncate(vector<vector<double>> mat, double index);
vector<double> MatVecProd_simple(vector<vector<double>>& mat , vector<double>& vec);
vector<double> SolveLinearDiag(vector<double> D, vector<double>z );
vector<double> forSub(vector<vector<double>> mat, vector<double>b);
vector<vector<double>> extractLowTri(vector<double> JM,vector<double> VM);
vector<vector<double>> extractMat(vector<double> JM,vector<double> VM);

//------------------------------------------------------------------------------------------------------------------------


// function that extracts lower triangle out of JM and VM arrays
vector<vector<double>> extractLowTri(vector<double> JM,vector<double> VM)
{
    int size = JM[0]-2;
    vector<vector<double>> low = CreateMatrix(size,size);
    int VMID = size;
    for(int i= 0; i < size ; i++)
    {
        low[i][i] = VM[i];

        if(JM[i+1]==JM[i])
        {
            i++;
            low[i][i] = VM[i];
        }

        double values = JM[i+1] - JM[i];   //number of values per row

        for(int j = 0; j< values; j++)
        {
            int colID = JM[VMID+1+j]-1; //-1 -> indexshift

            if (i>colID)
            {
            low[i][colID] = VM[VMID+1+j];
            }

            low[i][i] = VM[i];
        }
        VMID += values;

    }
    return low;
}

//------------------------------------------------------------------------------------------------------------------------
// function that extracts matrix out of JM and VM arrays
vector<vector<double>> extractMat(vector<double> JM,vector<double> VM)
{
    int size = JM[0]-2;
    vector<vector<double>> mat = CreateMatrix(size,size);
    int VMID = size;
    for(int i= 0; i < size ; i++)
    {
        mat[i][i] = VM[i];

        if(JM[i+1]==JM[i])
        {
            i++;
            mat[i][i] = VM[i];
        }

        double values = JM[i+1] - JM[i];   //number of values per row

        for(int j = 0; j< values; j++)
        {
            int colID = JM[VMID+1+j]-1; //-1 -> indexshift


            mat[i][colID] = VM[VMID+1+j];


            mat[i][i] = VM[i];
        }
        VMID += values;

    }
    return mat;
}

//------------------------------------------------------------------------------------------------------------------------

// solves a linear system where Matrix A is just a diagonal (stored in a vector)
//Ax=b => Mw=z => solve for w
vector<double> SolveLinearDiag(vector<double> D, vector<double>z )
{
    vector<double> w;
    for(int i = 0; i< D.size(); i++)
    {
        w.push_back(z[i] / D[i]);
    }
    return w;
}

//------------------------------------------------------------------------------------------------------------------------

// truncates the input !square!matrix up to the given index (indices start at 0)
vector<vector<double>> MatrixTruncate(vector<vector<double>> mat, double index)
{
    vector<vector<double>> new_mat;
    for(int i = 0; i <= index; i++)
    {
        vector<double> tmp = VecSlice(mat[i],0,index);
        new_mat.push_back(tmp);
    }
    return new_mat;
}

//------------------------------------------------------------------------------------------------------------------------

// solves Ax=b for A being a lower triangular matrix
vector<double> forSub(vector<vector<double>> mat, vector<double>b)
{
    int size = mat.size();
    vector<double> x(size);
    if(size != b.size())
    {
        cout << "Dimension of Matrix and Vector do not fit together! Programm terminates now!";
        exit(0);
    }

    for(int i = 0; i< size ; i++)
    {
        double sum = 0;
        for(int j = 0; j<i; j++)
        {
            sum += mat[i][j] * x[j];
        }
        x[i] = (b[i]-sum) / mat[i][i];
    }

    return x;
}

//------------------------------------------------------------------------------------------------------------------------

// solves Ax=b for A being a !strict upper triangular matrix!
vector<double> backSub(vector<vector<double>>& mat , vector<double>& b)
{
    int size = mat.size();
    vector<double> x(size,0);

    for(int i = size-1; i>= 0; i--)
    {
        double sum = 0;
        for(int j = i +1; j<size; j++)
        {
            sum += mat[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / mat[i][i];
    }
    return x;
}

//------------------------------------------------------------------------------------------------------------------------

// returns transposed of input matrix
vector<vector<double>> MatrixTrans(vector<vector<double>>& mat)
{
    double row = mat.size();
    double coloumn = mat[0].size();
    vector<vector<double>> trans = CreateMatrix(coloumn,row);
    for(int i = 0; i<row; i++)
    {
        for(int j = 0; j < coloumn; j++)
        {
            trans[j][i] = mat[i][j];
        }
    }
    return trans;
}

//------------------------------------------------------------------------------------------------------------------------

vector<double> VecScalarMult(vector<double> vec , double scalar)
{
    vector<double> result;
    for(int i=0;i<vec.size();i++)
        result.push_back(vec[i]*scalar);
    return result;
}

//------------------------------------------------------------------------------------------------------------------------

// performs Matrix Vector product for Matrices stored in "classic" format
// also works for nonsquare matrices
vector<double> MatVecProd_simple(vector<vector<double>>& mat , vector<double>& vec)
{
    if (mat[0].size() != vec.size()) // check if Dimensions fit
    {
        cout << "Dimension of Matrix and Vector do not fit together! Programm terminates now!";
        exit(0);
    }
    vector<double> result;
    for(int i = 0; i<mat.size(); i++)
    {
        result.push_back(DotProd(mat[i],vec));
    }
    return result;
}

//------------------------------------------------------------------------------------------------------------------------

// JM contains value pointers and coloumn indices
// VM contains diagonal elements and offdiagonal elements
vector<double> MatVecProd(bool& sym, std::vector <double>& JM , std::vector <double>& VM, std::vector<double>& x)
{
//first value of JM is n+2
int size = JM[0]-2; // size = n (A = nxn)
// create vectors
vector<double> A_D; //contains all diagonal elements of A
A_D = VecSlice(VM,0,size-1);

if (A_D.size() != x.size()) // check if Dimensions fit
{
    cout << "Dimension of Matrix and Vector do not fit together! Programm terminates now!";
    exit(0);
}

vector<double> A_off_D; //contains all offdiagonal elements of A. similar to (CSR) array V in tutorial 1
A_off_D = VecSlice(VM,size+1,VM.size()-1);

// solution vectors
vector<double> b(size,0);
vector<double> b_D; //solution vector for diagonal elements

//create seperate vector for the coloumn indices (offdiagonal values)
vector<double> J = VecSlice(JM,size+1,JM.size()-1);

// calculate diagonal solution part, this is similar for both cases (sym and not sym Matrices)
for(int i= 0; i<size; i++)
    {
        b_D.push_back(A_D[i]*x[i]);
    }

// The following calculation is done for both cases
for(int i = 0; i < size;i++)
    {
        double index_start = JM[i];
        //double index_start_x = ;         // in tutorial called J(i2)
        double index_end = JM[i+1];
        if(index_start == index_end)
            {
                b[i] = 0;
            }
        else
        {
            vector<double> VM_sliced = VecSlice(VM,index_start-1,index_end-2); //-1 because element 1 is at index 0
            vector<double> J_sliced = VecSlice(J,index_start-1-(size+1),index_end-2-(size+1));
            for(int s = 0; s<J_sliced.size(); s++)
                J_sliced[s]=J_sliced[s]-1; //index shift
            vector<double> x_sliced = PickSlice(x,J_sliced);

            b[i] = DotProd(VM_sliced,x_sliced);
        }
    }
b = VecAdd(b,b_D); // add diagonal and offidagonal results

vector<double> b_sym(size,0);
if(sym == 1)
    {
        for( int i = 0; i<size; i++)
        {
            double index_start = JM[i];
            double index_end = JM[i+1];
            if(index_start == index_end)
            {}
            else
            {
            vector<double> VM_sliced = VecSlice(VM,index_start-1,index_end-2); //-1 because element 1 is at index 0
            vector<double> J_sliced = VecSlice(J,index_start-1-(size+1),index_end-2-(size+1));
            for( int s = 0; s< J_sliced.size(); s++)
            {
                J_sliced[s]=J_sliced[s]-1; //index shift
                int rowid = J_sliced[s];
                b_sym[rowid] += VM_sliced[s]*x[rowid];
            }
            }
        }

    }
b = VecAdd(b,b_sym); // also add symmetric elements
return b;
}

//------------------------------------------------------------------------------------------------------------------------

//creates a matrix filled with 0s. Elements can be accessed by mat[i][j]
vector<vector<double>> CreateMatrix(double row, double coloumn)
{
    vector<vector<double>> mat(row, vector<double>(coloumn)); // creates matrix filled with 0s
    return mat;
}

//------------------------------------------------------------------------------------------------------------------------

// prints vector element by element
void PrintVec(vector<double>& vec)
{
        for (double element : vec)
    {
        cout << element << "  " << '\t';
    }
    cout << endl;
}

//------------------------------------------------------------------------------------------------------------------------

// prints a matrix row by row
void PrintMatrix(vector<vector<double>> mat)
{
    for(int i = 0; i <mat.size(); i++)
    {
        PrintVec(mat[i]);
    }
    cout << endl;
}

//------------------------------------------------------------------------------------------------------------------------

// adds two vectors together and returns resulting vector
vector<double> VecAdd(std::vector<double> vec1, std::vector<double> vec2)
    {
        std::vector<double> result;
        if(vec1.size() == vec2.size())
            {
                int size = vec1.size();
                for(int i = 0; i < size; i++)
                    {
                        result.push_back(vec1[i] + vec2[i]);
                    }
            }

        else
            {
                cout << "Vectors dont have the same size!" << endl;
            }
        return result ;
    }

//------------------------------------------------------------------------------------------------------------------------

// subtracts vector 2 from vector 1 and returns resulting vector
vector<double> VecSub(vector<double> vec1, vector<double> vec2)
    {
        vector<double> result;
        if(vec1.size() == vec2.size())
            {
                int size = vec1.size();
                for(int i = 0; i < size; i++)
                    {
                        result.push_back(vec1[i] - vec2[i]);
                    }
            }

        else
            {
                cout << "Vectors dont have the same size!" << endl;
            }
        return result ;
    }

//------------------------------------------------------------------------------------------------------------------------


// function that picks out certain elements of a vector, the indices of the elements need to be provided in a vector
vector<double> PickSlice(vector<double> vec,vector<double> indices)
    {
        vector<double> slice;

        for(int i = 0; i< indices.size(); i++)
        {
            double vec_index = indices[i];
            slice.push_back(vec[vec_index]);
        }
        return slice;
    }

//------------------------------------------------------------------------------------------------------------------------

//calculates dotproduct of two given vectors
double DotProd(std::vector<double> vec1, std::vector<double> vec2)
    {
        double dot_prod = 0;
        for(int i = 0; i < vec1.size(); i++)
            dot_prod += vec1[i] * vec2[i];
        return dot_prod;
    }

//------------------------------------------------------------------------------------------------------------------------

// function that slices a vector, including start point, end point and all elements inbetween
// also works when same index is inserted twice or when vector is not initialized or reserved
vector<double> VecSlice(vector<double>& vec,
                    int start, int end)
{
    // Starting and Ending iterators
    auto start_it = vec.begin() + start;
    auto end_it = vec.begin() + end + 1;

    // To store the sliced vector
    std::vector<double> result(end - start + 1);

    // Copy vector using copy function()
    copy(start_it, end_it, result.begin());

    // Return the final sliced vector
    return result;
}

//------------------------------------------------------------------------------------------------------------------------


// function that takes A, x_0 and b. Output is the residual
// Ax = b  =>  r = b - Ax_0
vector<double> residualCalc(bool& sym,vector<double>& JM, vector<double>& VM, vector<double>& x, vector<double>& b)
{
    vector <double> product = MatVecProd(sym,JM,VM,x);
    vector<double> res = VecSub(b,product);
    return res;
}

//------------------------------------------------------------------------------------------------------------------------

// calculates the euclidian norm of a vector
double VecNorm(vector<double> vec)
{
    double tmp = 0;
    for(int i = 0 ; i<vec.size(); i++)
        tmp += vec[i]*vec[i];//pow(vec[i],2.0); // raise element to power 2
    double result = sqrt(tmp);
    return result;
}

//------------------------------------------------------------------------------------------------------------------------

// function that divides each element of a vector by the input value
vector<double> VecElemDiv(vector<double> vec, double value)
{
    vector<double> result;
    for(int i=0; i<vec.size(); i++)
        result.push_back(vec[i] / value);
    return result;
}

//------------------------------------------------------------------------------------------------------------------------

// normalizes a vector to length 1
vector<double> VecNormalize(vector<double> vec)
{
    double norm = VecNorm(vec);
    vector<double> result = VecElemDiv(vec,norm);
    return result;
}

//------------------------------------------------------------------------------------------------------------------------
