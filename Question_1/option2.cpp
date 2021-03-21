#include <iostream>

double* AllocateVector(int noRows)
{
    double* vector;
    vector = new double [noRows];
    return vector;
}

void DeallocateVector(double* vector)
{
    delete[] vector;
}

void PrintVector(int noRows, double* vector)
{
    for (int i = 0; i < noRows; i++)
    {
        if (noRows == 1)
        {
            std::cout << "[ " << vector[i] << " ]";
        }
        else if (i == 0)
        {
            std::cout << "[ " << vector[i] << " ";
        }
        else if (i == noRows - 1) 
        {
            std::cout << vector[i] << " ]" << std::endl;
        }
        else
        {
            std::cout << vector[i] << " ";
        }   
    }
}

int main(int argc, char* argv[]){
    //  get user to input length and n points
    int n;
    double length;
    std::cout << "specify length and n --" << std::endl;
    std::cout << "length << ";
    std::cin >> length;
    std::cout << "n << ";
    std::cin >> n;

    //  calculate and store step size in h 
    double h = length/((double)n);

    /* allocate memory for diagonals and assign d = main diagonal, u = upper
    diagonal, l = lower diagonal     */
    double* d = AllocateVector(n);
    double* u = AllocateVector(n);
    double* l = AllocateVector(n);

    // for loop to generate coefficients
    for (int i = 0; i < n; i++)
    {
        d[i] = 4.0;
        if (i == 0)
        {
            l[i] = 0.0;
            u[i] = 2.0;
        } else if (i == n-1)
        {
            l[i] = 2.0;
            u[i] = 0.0;
        } else 
        {
            l[i] = 1.0;
            u[i] = 1.0;
        }
    }
    
    //  print vectors for upper diagonal, main diagonal and lower diagonal
    PrintVector(n, u);
    PrintVector(n, d);
    PrintVector(n, l);

    //  deallocate memory for upper diagonal, main diagonal and lower diagonal
    DeallocateVector(l);
    DeallocateVector(u);
    DeallocateVector(d);
}