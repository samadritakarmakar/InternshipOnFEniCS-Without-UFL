/*
 * Authored by Samadrita Karmakar June 2019
 * As of Now supports std::vector<Eigen::Matrix<double... and double
 * Each element of a vector reprensents a time step
 */

#ifndef WRITETOOCTAVEASCII_HPP
#define WRITETOOCTAVEASCII_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

void WriteToOctaveASCII(std::string FileName, std::vector<double>& Data);

template<int rowSize, int colSize>
void WriteToOctaveASCII(std::string FileName, std::vector<Eigen::Matrix<double, rowSize, colSize>>& Data)
{
    std::ofstream asciiFile;
    asciiFile.open(FileName.c_str());
    for (int row=0; row<rowSize; row++)
    {
        for (int Step=0; Step<Data.size(); Step++)
        {
            for (int col=0; col<colSize; col++)
            {
                asciiFile<<Data[Step](row, col)<<" ";
            }
            asciiFile<<"\t";
        }
        asciiFile<<"\n";
    }
}

void WriteToOctaveASCII(std::string FileName, std::vector<double>& Data)
{
    std::ofstream asciiFile;
    asciiFile.open(FileName.c_str());
    for (int Step=0; Step<Data.size(); Step++)
    {
        asciiFile<<Data[Step]<<"\t";
    }
}

#endif // WRITETOOCTAVEASCII_HPP
