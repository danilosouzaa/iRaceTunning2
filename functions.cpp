#include "functions.hpp"

int contCasasDecimais(double num)
{
    int i = 0;
    double m = 10;
    double aux;
    for (i = 1; i <= 1000; i++)
    {
        aux = num;
        aux *= m;
        if (ceil(aux - 1e-6) == floor(aux + 1e-6))
        {
            return i;
        }
        m *= 10;
    }
    return -1;
}

/*calculates the maximum common divisor for two integers.*/
int cutMaxDivisorCommonRec(int m, int n)
{

    int t = 0;
    m = m < 0 ? -m : m; /* abs(u) */
    n = n < 0 ? -n : n;
    if (m < n)
    {
        t = m;
        m = n;
        n = t;
    }

    if (n == 0)
        return m;
    else
    {
        int resto = m % n;
        try
        {
            if ((n == INT_MAX) || (n == INT_MIN))
            {
                return 1;
            }
            //printf("%d %d - ", n,resto);
            return cutMaxDivisorCommonRec(n, resto);
        }
        catch (int e)
        {
            return 1;
        }
    }
}

int cutMaxDivisorCommonVector(std::vector<double> coefs, int nElem)
{
    int n = (int)coefs[nElem];
    int mdc = 1;
    while (nElem > 0)
    {
        int m = (int)coefs[nElem - 1];
        try
        {
            mdc = cutMaxDivisorCommonRec(m, n);
        }
        catch (int e)
        {
            mdc = 1;
        }
        n = mdc;
        if (n == 1)
        {
            return 1;
        }
        nElem--;
    }
    int m = (int)coefs[nElem];
    try
    {
        mdc = cutMaxDivisorCommonRec(m, n);
    }
    catch (int e)
    {
        mdc = 1;
    }
    n = mdc;
    return n;
}

int heuristicConstraintsCoefFrac(std::vector<double> &temp, int szConstraints)
{

    int i, aux = 0;

    int *contCasas = (int *)malloc(sizeof(int) * (szConstraints + 1));
    int contFrac = 0;
    int maior = -1;
    double maiorCoef = temp[0];
    for (i = 0; i < szConstraints + 1; i++)
    {
        //temp[i] = input[i];
        if (temp[i] > maiorCoef)
        {
            maiorCoef = temp[i];
        }
        if (ceil(temp[i] - 1e-6) == floor(temp[i] + 1e-6))
        {
            contCasas[i] = 0;
        }
        else
        {
            contFrac++;
            contCasas[i] = contCasasDecimais(temp[i]);
        }
        if (contCasas[i] > maior)
        {
            maior = contCasas[i];
        }
    }
    //temp[aux] = constraintsOriginal->rightSide[idContraints];
    free(contCasas);

    // if (ceil(temp[szConstraints] - 1e-6) == floor(temp[szConstraints] + 1e-6))
    // {
    //   typeTemp[szConstraints] = 0;
    //   contCasas[szConstraints] = 0;
    // }
    // else
    // {
    //   contFrac++;
    //   typeTemp[szConstraints] = 1;
    //   contCasas[szConstraints] = contCasasDecimais(temp[szConstraints]);
    // }
    // if (contCasas[szConstraints] > maior)
    // {
    //   maior = contCasas[szConstraints];
    // }
    //printf("MAIORRR!!!: %d %f %d\n",maior, maiorCoef,(INT_MAX/20));
    if (((int)(pow(10, maior) * maiorCoef) <= (INT_MAX / 20)) && ((int)(pow(10, maior) * maiorCoef) >= -1 * (INT_MAX / 20)))
    {
        //input.clear();
        double m = pow(10, maior);
        for (i = 0; i < szConstraints + 1; i++)
        {
            temp[i] = (int)(temp[i] * m);
        }
        int mdc;
        try
        {
            mdc = cutMaxDivisorCommonVector(temp, szConstraints);
        }
        catch (int e)
        {
            mdc = 1;
        }
        // printf("mdc %d\n", mdc);

        //  printf("converteu 1: ");
        //  for (i = 0; i < szConstraints + 1; i++)
        //  {
        //    printf("%lf ", temp[i]);
        //  }
        // getchar();
        for (i = 0; i < szConstraints + 1; i++)
        {
            if (mdc > 1)
            {
                temp[i] /= mdc;
            }
            if (i < szConstraints)
            {
                temp[i] = floor(temp[i]);
            }
            else
            {
                temp[i] = ceil(temp[i]);
            }
            // input.push_back(temp[i]);
        }
        // printf("\nconverteu: ");
        // for (i = 0; i < szConstraints + 1; i++)
        // {
        //   printf("%lf ", temp[i]);
        // }

        return 1;
    }
    else
    {
        return 0;
    }
}

double valueViolation(constraintsReal *cCover, constraintsReal *constraintsSmall, int idCover, int ogConstraint, int precision)
{
    int i = 0, el;
    double lhs = 0, violation = 0;
    int aux = constraintsSmall->ElementsConstraints[ogConstraint];
    for (i = cCover->ElementsConstraints[idCover]; i < cCover->ElementsConstraints[idCover + 1]; i++)
    {
        el = constraintsSmall->Elements[aux];
        //lhs += ((double)cCover->Coefficients[i]) * ((double)constraintsSmall->xAsterisc[el]/(double)precision);
        lhs += cCover->Coefficients[i] * constraintsSmall->xAsterisc[el];
        //       printf("x* %d= %d coef = %d \n",el, constraintsSmall->xAsterisc[el], constraintsSmall->Coefficients[i]);
        aux++;
        //    }
    }
    //printf("lhs %f rhs %f\n",lhs, cCover->rightSide[idCover]);
    violation = lhs - (double)cCover->rightSide[idCover];
    if (violation + 1e-5 < 0.0)
    {
        violation = 0.0;
    }
    return violation;
}

int verifyRepeatedIncidency(int **matrizIncidencia, constraintsReal *originalConstraints, int posCover)
{
    int i, j, cont = 0, el, aux;
    int sz = originalConstraints->ElementsConstraints[posCover + 1] - originalConstraints->ElementsConstraints[posCover];
    for (i = 0; i < posCover; i++)
    {
        if ((sz == (originalConstraints->ElementsConstraints[i + 1] - originalConstraints->ElementsConstraints[i])) && (originalConstraints->rightSide[posCover] == originalConstraints->rightSide[i]))
        {
            for (j = originalConstraints->ElementsConstraints[posCover]; j < originalConstraints->ElementsConstraints[posCover + 1]; j++)
            {
                el = originalConstraints->Elements[j];
                if (matrizIncidencia[i][el] != -1)
                {
                    aux = matrizIncidencia[i][el];
                    if (originalConstraints->Coefficients[j] == originalConstraints->Coefficients[aux])
                    {
                        cont++;
                        //printf("cont %d sz %d\n", cont, sz);
                    }
                }
                else
                {
                    break;
                }
            }
        }
        if ((cont == sz) && (originalConstraints->rightSide[i] == originalConstraints->rightSide[posCover]))
        {
            return 0;
        }
        else
        {
            cont = 0;
        }
    }
    return 1;
}

constraintsReal *createCutsCover(constraintsReal *cutsCover, constraintsReal *constraintsOriginal, constraintsReal *constraintsSmall, int constraint, int nCuts, int precision)
{

    if (nCuts <= 0)
    {
        return constraintsOriginal;
    }
    long int i = 0, j = 0, cont = 0, contConstraints = 0, cTest = 0;
    ;
    int *idc_Cover = (int *)malloc(sizeof(int) * nCuts);
    for (i = 0; i < nCuts; i++)
    {
        // if (idc_Cover[i] == 0)
        // {
        //     continue;
        // }
        idc_Cover[i] = 1;
        double violation = valueViolation(cutsCover, constraintsSmall, i, constraint, precision);
        //printf("%f\n",violation);
        if (violation <= 0)
        {
            idc_Cover[i] = 0;
        }
        if (idc_Cover[i] == 1)
        {
            contConstraints++;
            cTest = 0;
            for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
            {
                if (cutsCover->Coefficients[j] != 0)
                {
                    cont++;
                    cTest++;
                }
            }
            if (cTest < 2)
            {
                cont -= cTest;
                contConstraints--;
                idc_Cover[i] = 0;
            }
        }
    }
    if (contConstraints == 0)
    {
        free(idc_Cover);
        return constraintsOriginal;
    }
    //printf("cont: %d contConstraints: %d original cont: %d original const:%d nCuts: %d \n", cont, contConstraints, constraintsOriginal->cont, constraintsOriginal->numberConstraints, nCuts);
    constraintsReal *outCutsNew = AllocStrConstraintsReal(constraintsOriginal->cont + cont, constraintsOriginal->numberConstraints + contConstraints, constraintsOriginal->numberVariables);
    outCutsNew->ElementsConstraints[0] = constraintsOriginal->ElementsConstraints[0];
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        outCutsNew->rightSide[i] = constraintsOriginal->rightSide[i];
        outCutsNew->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        outCutsNew->Elements[i] = constraintsOriginal->Elements[i];
        outCutsNew->Coefficients[i] = constraintsOriginal->Coefficients[i];
    }
    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        outCutsNew->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
        outCutsNew->ub[i] = constraintsOriginal->ub[i];
        outCutsNew->lb[i] = constraintsOriginal->lb[i];
        outCutsNew->typeVariables[i] = constraintsOriginal->typeVariables[i];
    }
    int aux = constraintsOriginal->numberConstraints;
    int c_aux = constraintsOriginal->cont;
    for (i = 0; i < nCuts; i++)
    {
        if (idc_Cover[i] == 1)
        {
            int c_XSolution = constraintsSmall->ElementsConstraints[constraint];
            outCutsNew->rightSide[aux] = cutsCover->rightSide[i]; // algo nessa linha
            for (j = cutsCover->ElementsConstraints[i]; j < cutsCover->ElementsConstraints[i + 1]; j++)
            {
                if (cutsCover->Coefficients[j] != 0)
                {
                    outCutsNew->Elements[c_aux] = constraintsSmall->Elements[c_XSolution]; //algo nessa linha
                    outCutsNew->Coefficients[c_aux] = cutsCover->Coefficients[j];
                    c_aux++;
                }
                c_XSolution++;
            }
            outCutsNew->ElementsConstraints[aux + 1] = c_aux;
            aux++;
        }
    }
    //printf("outCutNes: %d\n", outCutsNew->numberConstraints);
    int **matrizIncidencia = (int **)malloc(sizeof(int *) * outCutsNew->numberConstraints);
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        matrizIncidencia[i] = (int *)malloc(sizeof(int) * outCutsNew->numberVariables);
        for (j = 0; j < outCutsNew->numberVariables; j++)
        {
            matrizIncidencia[i][j] = -1;
        }
        for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
        {
            aux = outCutsNew->Elements[j];
            matrizIncidencia[i][aux] = j; //outCutsNew->ElementsConstraints[i];
        }
    }
    free(idc_Cover);
    int *validated = (int *)malloc(sizeof(int) * outCutsNew->numberConstraints);
    //validated[0] = 1;
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        if (i < constraintsSmall->numberConstraints)
        {
            validated[i] = 1;
        }
        else
        {
            validated[i] = verifyRepeatedIncidency(matrizIncidencia, outCutsNew, i);
        }
        // printf(" val. %d\t", validated[i]);
    }
    // getchar();
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        free(matrizIncidencia[i]);
    }
    free(matrizIncidencia);

    cont = 0, contConstraints = 0;

    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        if (validated[i] == 1)
        {
            contConstraints++;
            for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
            {
                cont++;
            }
        }
    }
    if (contConstraints == outCutsNew->numberConstraints)
    {
        // printf("t1\n");
        freeStrConstraintsReal(constraintsOriginal);
        free(validated);
        return outCutsNew;
    }
    constraintsReal *cutsNewNoRepetead = AllocStrConstraintsReal(cont, contConstraints, outCutsNew->numberVariables);
    cont = 0, contConstraints = 0;
    cutsNewNoRepetead->ElementsConstraints[0] = 0;
    for (i = 0; i < outCutsNew->numberConstraints; i++)
    {
        if (validated[i] == 1)
        {
            cutsNewNoRepetead->rightSide[contConstraints] = outCutsNew->rightSide[i];
            for (j = outCutsNew->ElementsConstraints[i]; j < outCutsNew->ElementsConstraints[i + 1]; j++)
            {
                cutsNewNoRepetead->Elements[cont] = outCutsNew->Elements[j];
                cutsNewNoRepetead->Coefficients[cont] = outCutsNew->Coefficients[j];
                cont++;
            }
            cutsNewNoRepetead->ElementsConstraints[contConstraints + 1] = cont;
            contConstraints++;
        }
    }
    for (i = 0; i < cutsNewNoRepetead->numberVariables; i++)
    {
        cutsNewNoRepetead->xAsterisc[i] = outCutsNew->xAsterisc[i];
    }

    free(validated);
    // printf("bug2");
    freeStrConstraintsReal(constraintsOriginal);
    freeStrConstraintsReal(outCutsNew);
    return cutsNewNoRepetead;
}

void quicksortCof(double *values, int *idc, int began, int end)
{
    int i, j;
    double pivo, aux;
    i = began;
    j = end - 1;
    pivo = values[(began + end) / 2];
    while (i <= j)
    {
        while (values[i] > pivo && i < end)
        {
            i++;
        }
        while (values[j] < pivo && j > began)
        {
            j--;
        }
        if (i <= j)
        {
            aux = values[i];
            values[i] = values[j];
            values[j] = aux;

            aux = idc[i];
            idc[i] = idc[j];
            idc[j] = aux;

            i++;
            j--;
        }
    }
    if (j > began)
        quicksortCof(values, idc, began, j + 1);
    if (i < end)
        quicksortCof(values, idc, i, end);
}

constraintsReal *AllocStrConstraintsReal(int cont, int nConstraints, int nVariables)
{
    constraintsReal *constraints = (constraintsReal *)malloc(sizeof(constraintsReal));
    constraints->Coefficients = (double *)malloc(sizeof(double) * cont);
    constraints->Elements = (int *)malloc(sizeof(int) * cont);
    constraints->ElementsConstraints = (int *)malloc(sizeof(int) * (nConstraints + 1));
    constraints->rightSide = (double *)malloc(sizeof(double) * nConstraints);
    constraints->xAsterisc = (double *)malloc(sizeof(double) * nVariables);
    constraints->typeVariables = (int *)malloc(sizeof(int) * nVariables);
    constraints->ub = (double *)malloc(sizeof(double) * nVariables);
    constraints->lb = (double *)malloc(sizeof(double) * nVariables);
    constraints->numberVariables = nVariables;
    constraints->numberConstraints = nConstraints;
    constraints->cont = cont;
    return constraints;
}

void showStructFull(constraintsReal *constraintsFull, std::vector<std::string> nameConstraints, std::vector<std::string> nameVariables)
{
    int i, j, el, flag;
    for (i = 0; i < constraintsFull->numberConstraints; i++)
    {
        flag = 0;
        std::cout << nameConstraints[i] << ": ";
        //printf("%s:\t ", nameConstraints[i]);
        for (j = constraintsFull->ElementsConstraints[i]; j < constraintsFull->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsFull->Elements[j];
            if (constraintsFull->Coefficients[j] > 0)
            {
                if (flag == 0)
                {
                    std::cout << constraintsFull->Coefficients[j] << " " << nameVariables[el] << "=" << constraintsFull->xAsterisc[el];
                    //printf("%.2f %s  = %lf ", constraintsFull->Coefficients[j], nameVariables[el], constraintsFull->xAsterisc[el]);
                    flag = 1;
                }
                else
                {
                    std::cout << " + " << constraintsFull->Coefficients[j] << " " << nameVariables[el] << "=" << constraintsFull->xAsterisc[el];
                    //printf("+ %.2f %s  = %lf ", constraintsFull->Coefficients[j], nameVariables[el], constraintsFull->xAsterisc[el]);
                }
            }
            else
            {
                flag = 1;
                std::cout << constraintsFull->Coefficients[j] << " " << nameVariables[el] << "=" << constraintsFull->xAsterisc[el];
                //printf("%.2f %s = %lf ", constraintsFull->Coefficients[j], nameVariables[el], constraintsFull->xAsterisc[el]);
            }
        }
        std::cout << " <= " << constraintsFull->rightSide[i] << std::endl;
        //printf("<= %.2f\n", constraintsFull->rightSide[i]);
    }
}

constraintsReal *prepareConstraintsOfLifting(constraintsReal *constraintsInitial, int precision)
{
    int *binaryConstraints = (int *)malloc(sizeof(int) * constraintsInitial->numberConstraints);
    int i, j, qntBin, qntNBin, el, auxc = 0, auxc2 = 0;
    for (i = 0; i < constraintsInitial->numberConstraints; i++)
    {
        qntBin = 0;
        qntNBin = 0;
        for (j = constraintsInitial->ElementsConstraints[i]; j < constraintsInitial->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsInitial->Elements[j];
            //printf("%d ", constraintsInitial->typeVariables[el]);
            if (constraintsInitial->typeVariables[el] == 0)
            {
                qntBin++;
            }
            else
            {
                qntNBin++;
            }
        }
        //getchar();
        if (qntBin <= 1)
        {
            binaryConstraints[i] = 0; // restrição que não vale a pena transformar
        }
        else if (qntNBin == 0)
        {
            binaryConstraints[i] = 1; // restrição com todos os elementos binários
            auxc2++;
        }
        else
        {
            binaryConstraints[i] = 2; // restrição vale a pena transformar em binária
            auxc++;
        }
    }
    int contConstraints = 0, newCont = 0, contAux;
    //printf("Total: %d total binario: %d total q pode ser aproveitado: %d\n", constraintsInitial->numberConstraints, auxc2, auxc);
    int *constraintsAvaliable = (int *)malloc(sizeof(int) * constraintsInitial->numberConstraints);
    for (i = 0; i < constraintsInitial->numberConstraints; i++)
    {
        constraintsAvaliable[i] = 0;
        if (binaryConstraints[i] != 0)
        {
            contAux = 0;
            for (j = constraintsInitial->ElementsConstraints[i]; j < constraintsInitial->ElementsConstraints[i + 1]; j++)
            {
                el = constraintsInitial->Elements[j];
                if (constraintsInitial->typeVariables[el] != 0)
                {
                    if (((constraintsInitial->Coefficients[j] + 1e-5 < 0) && (constraintsInitial->ub[el] == DBL_MAX)) || ((constraintsInitial->Coefficients[j] - 1e-5 > 0) && (constraintsInitial->lb[el] == -DBL_MAX)))
                    {
                        contAux = 0;
                        break;
                    }
                }
                else
                {
                    contAux++;
                }
            }
            if (contAux >= 2)
            {
                newCont += contAux;
                contConstraints++;
                constraintsAvaliable[i] = 1;
            }
        }
    }
    constraintsReal *constraintsBinary = AllocStrConstraintsReal(newCont, contConstraints, constraintsInitial->numberVariables);
    newCont = 0;
    contConstraints = 0;
    constraintsBinary->ElementsConstraints[0] = 0;
    int rhsTemp = 0;
    for (i = 0; i < constraintsInitial->numberConstraints; i++)
    {
        if (constraintsAvaliable[i] == 1)
        {
            rhsTemp = 0;
            for (j = constraintsInitial->ElementsConstraints[i]; j < constraintsInitial->ElementsConstraints[i + 1]; j++)
            {
                el = constraintsInitial->Elements[j];
                if (constraintsInitial->typeVariables[el] == 0)
                {
                    constraintsBinary->Coefficients[newCont] = constraintsInitial->Coefficients[j];
                    constraintsBinary->Elements[newCont] = el;
                    newCont++;
                }
                else
                {
                    if (constraintsInitial->Coefficients[j] < 0)
                    {
                        if ((constraintsInitial->ub[el] >= 0) && (constraintsInitial->lb[el] >= 0))
                        {
                            rhsTemp += (constraintsInitial->Coefficients[j] * -1) * constraintsInitial->ub[el];
                        }
                        else if ((constraintsInitial->ub[el] <= 0) && (constraintsInitial->lb[el] <= 0))
                        {
                            rhsTemp += 0;
                        }
                        else
                        {
                            rhsTemp += (constraintsInitial->Coefficients[j] * -1) * constraintsInitial->ub[el];
                        }
                    }
                    else
                    {
                        if ((constraintsInitial->ub[el] >= 0) && (constraintsInitial->lb[el] >= 0))
                        {
                            rhsTemp += 0;
                        }
                        else if ((constraintsInitial->ub[el] <= 0) && (constraintsInitial->lb[el] <= 0))
                        {
                            rhsTemp += constraintsInitial->Coefficients[j] * (constraintsInitial->lb[el] * (-1));
                        }
                        else
                        {
                            rhsTemp += constraintsInitial->Coefficients[j] * (constraintsInitial->lb[el] * -1);
                        }
                    }
                }
            }
            constraintsBinary->ElementsConstraints[contConstraints + 1] = newCont;
            constraintsBinary->rightSide[contConstraints] = constraintsInitial->rightSide[i] + rhsTemp;
            contConstraints++;
        }
    }
    for (i = 0; i < constraintsInitial->numberVariables; i++)
    {
        //constraintsBinary->xAsterisc[i] = (int)(constraintsInitial->xAsterisc[i]*precision);
        constraintsBinary->xAsterisc[i] = constraintsInitial->xAsterisc[i];
        //printf("%lf \n", constraintsBinary->xAsterisc[i]);
        constraintsBinary->ub[i] = constraintsInitial->ub[i];
        constraintsBinary->lb[i] = constraintsInitial->lb[i];
        constraintsBinary->typeVariables[i] = constraintsInitial->typeVariables[i];
    }
    free(constraintsAvaliable);
    free(binaryConstraints);
    // printf("Total: %d\n", constraintsBinary->numberConstraints);
    return constraintsBinary;
}

void freeStrConstraintsReal(constraintsReal *cut)
{
    free(cut->Coefficients);
    free(cut->Elements);
    free(cut->ElementsConstraints);
    free(cut->rightSide);
    free(cut->xAsterisc);
    free(cut->lb);
    free(cut->ub);
    free(cut->typeVariables);
    free(cut);
}

constraintsReal *removeNegativeCoefficientsAndSort(constraintsReal *constraintsOriginal, int *convertVector)
{
    int i, j;
    int qntX = constraintsOriginal->numberVariables;
    int qntNegative = 0;
    double rhs = 0;
    int el = 0;
    for (i = 0; i < constraintsOriginal->cont; i++)
    {
        if (constraintsOriginal->Coefficients[i] < 0)
        {
            qntNegative++;
        }
    }
    qntNegative += constraintsOriginal->numberVariables;
    constraintsReal *newConstraints = AllocStrConstraintsReal(constraintsOriginal->cont, constraintsOriginal->numberConstraints, qntNegative);

    for (i = 0; i < constraintsOriginal->numberVariables; i++)
    {
        newConstraints->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
    }

    newConstraints->ElementsConstraints[0] = 0;
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        rhs = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            if (constraintsOriginal->Coefficients[j] < 0)
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j] * (-1);
                rhs += newConstraints->Coefficients[j];

                newConstraints->Elements[j] = qntX;
                el = constraintsOriginal->Elements[j];
                //printf("%lf\n", constraintsOriginal->xAsterisc[el]);
                newConstraints->xAsterisc[qntX] = 1 - constraintsOriginal->xAsterisc[el];
                convertVector[qntX - constraintsOriginal->numberVariables] = constraintsOriginal->Elements[j];
                qntX++;
            }
            else
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = constraintsOriginal->Elements[j];
            }
        }
        newConstraints->rightSide[i] = rhs;
        newConstraints->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
    }

    SortByCoefficients(newConstraints);
    freeStrConstraintsReal(constraintsOriginal);
    return newConstraints;
}

void SortByCoefficients(constraintsReal *h_cut)
{
    int i = 0;
    for (i = 0; i < h_cut->numberConstraints; i++)
    {
        quicksortCof(h_cut->Coefficients, h_cut->Elements, h_cut->ElementsConstraints[i], h_cut->ElementsConstraints[i + 1]);
    }
}

void quicksortDouble(double *values, int began, int end)
{
    int i, j;
    double pivo, aux;
    i = began;
    j = end - 1;
    pivo = values[(began + end) / 2];
    while (i <= j)
    {
        while (values[i] > pivo && i < end)
        {
            i++;
        }
        while (values[j] < pivo && j > began)
        {
            j--;
        }
        if (i <= j)
        {
            aux = values[i];
            values[i] = values[j];
            values[j] = aux;
            i++;
            j--;
        }
    }
    if (j > began)
        quicksortDouble(values, began, j + 1);
    if (i < end)
        quicksortDouble(values, i, end);
}

int *LCIBallas(int *coverSolution, constraintsReal *constraintsCover, int constraint)
{
    int i, w, szCover = 0;
    int szConstraints = constraintsCover->ElementsConstraints[constraint + 1] - constraintsCover->ElementsConstraints[constraint];
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        szCover += coverSolution[i - constraintsCover->ElementsConstraints[constraint]];
    }
    double *S_barra = (double *)malloc(sizeof(double) * (szCover + 1));
    S_barra[0] = 0;
    w = 1;
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        if (coverSolution[i - constraintsCover->ElementsConstraints[constraint]] == 1)
        {
            S_barra[w] = S_barra[w - 1] + constraintsCover->Coefficients[i];
            w++;
        }
    }

    int ini = 0;
    int *cutCoverLifted = (int *)malloc(sizeof(int) * (szConstraints + 1));
    for (w = constraintsCover->ElementsConstraints[constraint]; w < constraintsCover->ElementsConstraints[constraint + 1]; w++)
    {

        //cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
        if (coverSolution[w - constraintsCover->ElementsConstraints[constraint]] == 1)
        {
            cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = 1;
        }
        else
        {
            int flag = 0;
            for (ini = 0; ini < szCover; ini++)
            {
                if ((constraintsCover->Coefficients[w] > S_barra[ini]) && (constraintsCover->Coefficients[w] <= S_barra[ini + 1]))
                {
                    cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = ini;
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = constraintsCover->Coefficients[w];
            }
        }
    }

    cutCoverLifted[szConstraints] = szCover - 1;
    free(S_barra);
    return cutCoverLifted;
}

int *LCIAdam(int *coverSolution, constraintsReal *constraintsCover, int constraint)
{
    double fillBag = 0;
    int i, szCover = 0;
    int szConstraints = constraintsCover->ElementsConstraints[constraint + 1] - constraintsCover->ElementsConstraints[constraint];
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        if (constraintsCover->Coefficients[i] > constraintsCover->rightSide[constraint])
        {
            coverSolution[i - constraintsCover->ElementsConstraints[constraint]] = 0;
        }
        fillBag += (double)coverSolution[i - constraintsCover->ElementsConstraints[constraint]] * constraintsCover->Coefficients[i];
        szCover += coverSolution[i - constraintsCover->ElementsConstraints[constraint]];
    }
    if ((szCover <= 1) || (fillBag <= (double)constraintsCover->rightSide[constraint]))
    {
        return NULL;
    }
    int *n_coef = (int *)malloc(sizeof(int) * szCover);
    int *n_el = (int *)malloc(sizeof(int) * szCover);
    int caux = 0, qnt = 0;
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        if (coverSolution[caux] == 1)
        {
            n_coef[qnt] = (int)constraintsCover->Coefficients[i];
            n_el[qnt] = i;
            qnt++;
        }
        caux++;
    }
    double b = (double)constraintsCover->rightSide[constraint];
    double delta = 0;
    double phi = ((double)fillBag - (double)constraintsCover->rightSide[constraint]);
    int k = 1, w;
    double a_barra = (double)n_coef[0];
    for (w = 1; w < qnt; w++)
    {

        delta = a_barra - (double)n_coef[w];
        if (((double)k * delta) < phi)
        {
            a_barra = (double)n_coef[w];
            phi = phi - ((double)k * delta);
        }
        else
        {
            a_barra = a_barra - (phi / (double)k);
            phi = 0;
            break;
        }
        k++;
    }
    if (phi - 1e-5 > 0)
    {
        a_barra = (double)b / (double)szCover;
    }
    // printf("a_barra: %f\n", a_barra);
    int *c_menus = (int *)malloc(sizeof(int) * szCover);
    int *c_mais = (int *)malloc(sizeof(int) * szCover);
    double *S_barra = (double *)malloc(sizeof(double) * (szCover + 1));
    double *aux = (double *)malloc(sizeof(double) * szCover);
    int id1 = 0, id2 = 0;

    for (w = 0; w < qnt; w++)
    {
        if ((double)n_coef[w] - 1e-5 <= a_barra)
        {
            c_menus[id1] = w;
            id1++;
        }
        else
        {
            c_mais[id2] = w;
            id2++;
        }
    }
    for (w = 0; w < qnt; w++)
    {
        if ((double)n_coef[w] - 1e-5 < a_barra)
        {
            aux[w] = n_coef[w];
        }
        else
        {
            aux[w] = a_barra;
        }
    }
    quicksortDouble(aux, 0, qnt);
    S_barra[0] = 0;
    for (w = 1; w < qnt; w++)
    {
        S_barra[w] = S_barra[w - 1] + aux[w - 1];
    }
    S_barra[qnt] = constraintsCover->rightSide[constraint];

    int ini = 0;
    int *cutCoverLifted = (int *)malloc(sizeof(int) * (szConstraints + 1));
    for (w = constraintsCover->ElementsConstraints[constraint]; w < constraintsCover->ElementsConstraints[constraint + 1]; w++)
    {
        int flag = 0;
        //cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
        for (ini = 0; ini < qnt; ini++)
        {
            // printf("%lf %lf\n", constraintsCover->Coefficients[w],S_barra[ini]);
            if ((constraintsCover->Coefficients[w] > S_barra[ini]) && (constraintsCover->Coefficients[w] <= S_barra[ini + 1]))
            {
                cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = ini;
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            //printf("caiu aqui: %d s_ %f %f %f \n", qnt, S_barra[qnt], constraintsCover->Coefficients[w], constraintsCover->rightSide[constraint]);
            //getchar();
            cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = 0; //constraintsCover->Coefficients[w];
        }
    }
    int el;
    for (w = 0; w < id1; w++)
    {
        el = n_el[c_menus[w]];
        cutCoverLifted[el - constraintsCover->ElementsConstraints[constraint]] = 1;
    }
    cutCoverLifted[szConstraints] = szCover - 1;
    // printf("\nLifting:");
    // for(i=0;i<szConstraints+1;i++){
    //   printf("%d ", cutCoverLifted[i]);
    // }
    // printf("\n");
    free(c_menus);
    free(c_mais);
    free(S_barra);
    free(aux);
    free(n_coef);
    free(n_el);
    return cutCoverLifted;
}

int verifyViolationGreedy(int *solutionFinal, constraintsReal *constraitsUsed, int constraint, int precision)
{
    double lhs = 0.0;
    int i, el;
    int sz = constraitsUsed->ElementsConstraints[constraint + 1] - constraitsUsed->ElementsConstraints[constraint];
    for (i = constraitsUsed->ElementsConstraints[constraint]; i < constraitsUsed->ElementsConstraints[constraint + 1]; i++)
    {
        el = constraitsUsed->Elements[i];
        double aux = (double)solutionFinal[i - constraitsUsed->ElementsConstraints[constraint]]; //aqui era int ver o que acontece??
        aux *= constraitsUsed->xAsterisc[el];
        lhs += aux;
    }
    //lhs = lhs/(double)precision;;
    if (lhs + 1e-5 > solutionFinal[sz])
    {
        //printf("%f %d\n", lhs, solutionFinal[sz]);
        return 1;
    }
    return 0;
}

double avaliaSolution(constraintsReal *knapsackConstraints, int *solution, int constraint, int precision, int typeLift)
{
    int i, aux = 0;
    double lhs = 0.0;
    for (i = knapsackConstraints->ElementsConstraints[constraint]; i < knapsackConstraints->ElementsConstraints[constraint + 1]; i++)
    {
        lhs += (double)solution[aux] * knapsackConstraints->Coefficients[i];
        aux++;
    }
    if (lhs <= knapsackConstraints->rightSide[constraint])
    {
        return -99999;
    }
    int *solutionLCI;
    solutionLCI = LCIAdam(solution, knapsackConstraints, constraint);
    //solutionLCI = LCIBallas(solution, knapsackConstraints, constraint);
    if (solutionLCI == NULL)
    {
        return -99999;
    }
    double fo = 0;
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    for (i = knapsackConstraints->ElementsConstraints[constraint]; i < knapsackConstraints->ElementsConstraints[constraint + 1]; i++)
    {
        int el = knapsackConstraints->Elements[i];
        int p = i - knapsackConstraints->ElementsConstraints[constraint];
        fo += (double)solutionLCI[p] * knapsackConstraints->xAsterisc[el];
    }
    fo = fo - solutionLCI[sz];
    //fo /= precision;
    free(solutionLCI);
    return fo;
}

int *greedyInitialSolutionSA(constraintsReal *knapsackConstraints, int constraint, int precision, int typeLift, int typeSolutionInitial)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *solutionCover = (int *)malloc(sizeof(int) * sz);
    double *xTemp = (double *)malloc(sizeof(double) * sz);
    int *idc = (int *)malloc(sizeof(int) * sz);
    int i, el, aux = 0;
    for (i = knapsackConstraints->ElementsConstraints[constraint]; i < knapsackConstraints->ElementsConstraints[constraint + 1]; i++)
    {
        el = knapsackConstraints->Elements[i];
        if (typeSolutionInitial == 0)
        {
            //printf("With 1\n");
            //xTemp[aux] = (double)knapsackConstraints->xAsterisc[el] / (double)precision;
            xTemp[aux] = (double)knapsackConstraints->xAsterisc[el];
        }
        else
        {
            //printf("With 2\n");
            //xTemp[aux] = ((double)knapsackConstraints->xAsterisc[el] / (double)precision) / knapsackConstraints->Coefficients[i];
            xTemp[aux] = ((double)knapsackConstraints->xAsterisc[el] / knapsackConstraints->Coefficients[i]);
        }
        idc[aux] = i;
        solutionCover[aux] = 0;
        aux++;
    }
    //getchar();
    quicksortCof(xTemp, idc, 0, sz);
    int lhs = 0, posRef = 0;
    for (i = 0; i < sz; i++)
    {
        el = idc[i];
        lhs += knapsackConstraints->Coefficients[el];
        if (lhs > knapsackConstraints->rightSide[constraint])
        {
            lhs -= knapsackConstraints->Coefficients[el];
            posRef = i;
            break;
        }
        else
        {
            solutionCover[el - knapsackConstraints->ElementsConstraints[constraint]] = 1;
        }
    }
    int k = 0;
    for (k = posRef; k < sz; k++)
    {
        int *solutionFinal;
        el = idc[k];
        solutionCover[el - knapsackConstraints->ElementsConstraints[constraint]] = 1;
        lhs += knapsackConstraints->Coefficients[el];
        if (lhs <= knapsackConstraints->rightSide[constraint])
        {
            continue;
        }
        // for(i=0;i<sz;i++){
        //     printf("antes %d ", solutionCover[i]);
        // }
        //solutionFinal = LCIBallas(solutionCover, knapsackConstraints, constraint);
        solutionFinal = LCIAdam(solutionCover, knapsackConstraints, constraint);
        if (solutionFinal == NULL)
        {
            continue;
        }
        // printf("\n3\n");
        if (verifyViolationGreedy(solutionFinal, knapsackConstraints, constraint, precision) == 1)
        {
            free(xTemp);
            free(idc);
            free(solutionFinal);
            // for(i=0;i<sz;i++){
            //   printf("depois %d ", solutionCover[i]);
            // }
            return solutionCover;
        }
        free(solutionFinal);
    }
    free(xTemp);
    free(idc);
    free(solutionCover);
    return NULL;
}

void copyAndVerifyPoolSolution(int *solution, int sz, int *poolSolution, int *numberSolutionAtual, double violation, double *violationFinal)
{
    int i, j, atual;
    int flag = *numberSolutionAtual;
    atual = flag;
    for (i = 0; i < *numberSolutionAtual; i++)
    {
        for (j = 0; j < sz; j++)
        {
            if (solution[j] != poolSolution[j + i * sz])
            {
                flag--;
                break;
            }
        }
    }
    if (flag == 0)
    {
        *violationFinal += violation;
        for (j = 0; j < sz; j++)
        {
            poolSolution[j + atual * sz] = solution[j];
        }
        (*numberSolutionAtual)++;
    }
}

double fRand(double fMin, double fMax)
{
    struct timeval time;
    fflush(stdin);
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void shuffleVectorInt(int *vec, int sz)
{
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int j;
    double *num_temp = (double *)malloc(sizeof(double) * sz); //free
    for (j = 0; j < sz; j++)
    {
        num_temp[j] = (double)(rand() % RAND_MAX);
    }
    quicksortCof(num_temp, vec, 0, sz);
    free(num_temp);
}

//insert one remove one
int *neighborhoodOne(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }

    if (szNActive > 0)
    {
        if (szNActive > 1)
            shuffleVectorInt(vNActived, szNActive);
        newSolution[vNActived[0]] = 1;
    }
    if (szActive > 0)
    {
        if (szActive > 1)
            shuffleVectorInt(vActived, szActive);
        newSolution[vActived[0]] = 0;
    }
    free(vNActived);
    free(vActived);
    return newSolution;
}

//insert 2 remove 1
int *neighbohoodTwo(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }

    if (szNActive > 0)
    {
        if (szNActive > 1)
        {
            shuffleVectorInt(vNActived, szNActive);
            newSolution[vNActived[1]] = 1;
        }
        newSolution[vNActived[0]] = 1;
    }
    if (szActive > 0)
    {
        if (szActive > 1)
            shuffleVectorInt(vActived, szActive);
        newSolution[vActived[0]] = 0;
    }
    free(vNActived);
    free(vActived);
    return newSolution;
}

//Insert 1 remove 2
int *neighbohoodThree(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }

    //shuffleVectorInt(vActived, szActive);
    if (szNActive > 0)
    {
        if (szNActive > 1)
            shuffleVectorInt(vNActived, szNActive);
        newSolution[vNActived[0]] = 1;
    }

    if (szActive > 0)
    {
        if (szActive > 1)
        {
            shuffleVectorInt(vActived, szActive);
            newSolution[vActived[1]] = 0;
        }
        newSolution[vActived[0]] = 0;
    }

    free(vNActived);
    free(vActived);
    return newSolution;
}
//Insert 1
int *neighbohoodFour(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    //int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 0)
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    if (szNActive > 0)
    {
        if (szNActive > 1)
            shuffleVectorInt(vNActived, szNActive);
        newSolution[vNActived[0]] = 1;
    }
    free(vNActived);
    return newSolution;
}

//Remove 1
int *neighbohoodFive(constraintsReal *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    //int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
    }
    if (szActive > 0)
    {
        if (szActive > 1)
            shuffleVectorInt(vActived, szActive);
        newSolution[vActived[0]] = 0;
    }
    free(vActived);
    return newSolution;
}

constraintsReal *copyStrConstraintsReal(constraintsReal *originalRows)
{
    constraintsReal *destRows = AllocStrConstraintsReal(originalRows->cont, originalRows->numberConstraints, originalRows->numberVariables);
    int i;
    destRows->ElementsConstraints[0] = originalRows->ElementsConstraints[0];
    for (i = 0; i < originalRows->numberConstraints; i++)
    {
        destRows->rightSide[i] = originalRows->rightSide[i];
        destRows->ElementsConstraints[i + 1] = originalRows->ElementsConstraints[i + 1];
    }
    for (i = 0; i < originalRows->cont; i++)
    {
        destRows->Elements[i] = originalRows->Elements[i];
        destRows->Coefficients[i] = originalRows->Coefficients[i];
    }
    for (i = 0; i < originalRows->numberVariables; i++)
    {
        destRows->xAsterisc[i] = originalRows->xAsterisc[i];
        destRows->ub[i] = originalRows->ub[i];
        destRows->lb[i] = originalRows->lb[i];
        destRows->typeVariables[i] = originalRows->typeVariables[i];
    }
    return destRows;
}

constraintsReal *runCCwithSA(constraintsReal *binaryRows, int precision, int szPoolCutsMax, float TempInitial, float FactorResf, int typeLift, int minimal, double v1, double v2, double v3, double v4, double *violationFinal, int typeSolutionInitial, double timeLeft, int fm, double timeConstraints)
{

    //int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull); //verifica se as restrições podem ser transformadas em restrições da mochila
    //cutSmall *knapsackConstraints = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    //cutCover *coverCuts = CopyCutToCover(knapsackConstraints);
    int i, j, k;
    double randomico = 0.0;
    double delta_max = 0;
    double startT = omp_get_wtime();
    double currentTime = omp_get_wtime() - startT;
    //printf("antes: %d\n", binaryRows->numberConstraints);
    constraintsReal *copyRowsBinary = copyStrConstraintsReal(binaryRows);
    //printf("Depois: %d\n", binaryRows->numberConstraints);
    for (i = 0; i < copyRowsBinary->numberConstraints; i++)
    {
        double lhs = 0;
        int sz = copyRowsBinary->ElementsConstraints[i + 1] - copyRowsBinary->ElementsConstraints[i];
        int numberCuts = 0;
        double fo_best = 0;
        if ((copyRowsBinary->rightSide[i] <= 0))
        {
            //printf("Saiu por aqui: %d\n", i);
            continue;
        }

        int *bestSolution = greedyInitialSolutionSA(copyRowsBinary, i, precision, typeLift, typeSolutionInitial);
        if (bestSolution == NULL)
        {

            continue;
        }
        int *poolSolution = (int *)malloc(sizeof(int) * szPoolCutsMax * sz);
        // printf("\n2\n");
        // for(j=0;j<sz;j++){
        //   printf("durante: %d ", bestSolution[j]);
        // }
        fo_best = avaliaSolution(copyRowsBinary, bestSolution, i, precision, typeLift);
        // printf("violation Initial %lf\n", fo_best);
        if (fo_best - 1e-6 > 0)
        {
            copyAndVerifyPoolSolution(bestSolution, sz, poolSolution, &numberCuts, fo_best, violationFinal);
        }
        //printf("fo: %f\n", fo_asterisc);
        srand(time(NULL));
        double T = TempInitial;
        int *solutionVizinha;
        int semMelhora = sz * fm;
        int auxMelhora = 0;
        double f_asterisc = fo_best;
        double timeLeftConstraints = omp_get_wtime();
        //printf("tamanho: %d\n",sz);
        do
        {
            //break;//teste
            currentTime = omp_get_wtime() - startT;
            //printf("tempo: %lf %lf %d \n", timeLeft, currentTime, sz);
            if ((currentTime > timeLeft) || ((omp_get_wtime() - timeLeftConstraints) > timeConstraints))
            {

                break;
            }

            randomico = fRand(0.0, 1.0);

            if (randomico < v1)
            {

                solutionVizinha = neighborhoodOne(copyRowsBinary, precision, i, bestSolution);
            }
            else if (randomico < v1 + v2)
            {
                //printf("vz2");
                solutionVizinha = neighbohoodTwo(copyRowsBinary, precision, i, bestSolution);
            }
            else if (randomico < v1 + v2 + v3)
            {
                //printf("vz3");
                solutionVizinha = neighbohoodThree(copyRowsBinary, precision, i, bestSolution);
            }
            else if (randomico < v1 + v2 + v3 + v4)
            {
                solutionVizinha = neighbohoodFour(copyRowsBinary, precision, i, bestSolution);
            }
            else
            {
                solutionVizinha = neighbohoodFive(copyRowsBinary, precision, i, bestSolution);
            }
            double fo_vizinho = avaliaSolution(copyRowsBinary, solutionVizinha, i, precision, typeLift);
            // printf("violation vizinho: %f\n",fo_vizinho);
            if (fo_vizinho - 1e-6 > 0)
            {
                copyAndVerifyPoolSolution(solutionVizinha, sz, poolSolution, &numberCuts, fo_vizinho, violationFinal);
            }
            double delta = fo_vizinho - fo_best;
            if (delta > 0)
            {
                fo_best = fo_vizinho;
                for (j = 0; j < sz; j++)
                {
                    bestSolution[j] = solutionVizinha[j];
                }
                if (delta > delta_max)
                {
                    delta_max = delta;
                }
                if (fo_vizinho > f_asterisc)
                {
                    f_asterisc = fo_vizinho;
                    auxMelhora = 0;
                }
                else
                {
                    auxMelhora++;
                }
            }
            else
            {
                auxMelhora++;
                double x = fRand(0.0, 1);
                double v = exp(delta / T);
                if (x < v)
                {
                    fo_best = fo_vizinho;
                    for (j = 0; j < sz; j++)
                    {
                        bestSolution[j] = solutionVizinha[j];
                    }
                }
            }
            //printf("fo Vizinha:%f\n", fo_vizinho);
            T = T * FactorResf;
            free(solutionVizinha);

        } while ((T - 1e-5 > 0) && (auxMelhora <= semMelhora));
        //free(solutionVizinha);
        //*violationFinal += fo_asterisc;
        free(bestSolution);
        //printf("Delta Max: %f\n", delta_max);
        int c_AuxSolution = 0;
        int c_XSolution = 0;
        // printf("NumberCuts: %d\n", numberCuts);
        constraintsReal *cutsCoverSolution = AllocStrConstraintsReal(sz * numberCuts, numberCuts, copyRowsBinary->numberVariables);
        cutsCoverSolution->ElementsConstraints[0] = 0;
        int itePool;
        for (itePool = 0; itePool < numberCuts; itePool++)
        {
            //qnt = 0;
            int caux = 0;
            lhs = 0;
            int *solutionTemp = (int *)malloc(sizeof(int) * sz);
            for (k = copyRowsBinary->ElementsConstraints[i]; k < copyRowsBinary->ElementsConstraints[i + 1]; k++)
            {
                //  qnt += poolSolution[caux + itePool * szConstraint];

                solutionTemp[k - copyRowsBinary->ElementsConstraints[i]] = poolSolution[(k - copyRowsBinary->ElementsConstraints[i]) + itePool * sz];
                lhs += copyRowsBinary->Coefficients[k] * poolSolution[caux + itePool * sz];
                caux++;
            }
            if (lhs <= copyRowsBinary->rightSide[i])
            {
                free(solutionTemp);
                continue;
            }

            int *cutCoverLifted = LCIAdam(solutionTemp, copyRowsBinary, i);
            // int *cutCoverLifted = LCIBallas(solutionTemp, copyRowsBinary, i);
            //  printf("Restrição de origem:%d\n",i);
            //  for(k = copyRowsBinary->ElementsConstraints[i];k<copyRowsBinary->ElementsConstraints[i+1];k++){
            //    int el  =copyRowsBinary->Elements[k];
            //    printf("%lf x%d = %lf  + ", copyRowsBinary->Coefficients[k],el, copyRowsBinary->xAsterisc[el]);
            //  }
            //  printf("<=%lf\n",copyRowsBinary->rightSide[i]);
            for (k = 0; k < sz; k++)
            {
                cutsCoverSolution->Coefficients[c_XSolution] = (double)cutCoverLifted[k];
                int el = copyRowsBinary->ElementsConstraints[i] + k;
                cutsCoverSolution->Elements[c_XSolution] = copyRowsBinary->Elements[el];
                // printf("%d x%d = %d  + ", cutCoverLifted[k], copyRowsBinary->Elements[el],solutionTemp[k]);
                c_XSolution++;
            }
            // printf("<= %d\n", cutCoverLifted[sz]);
            cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;
            cutsCoverSolution->rightSide[c_AuxSolution] = (double)cutCoverLifted[sz];
            c_AuxSolution++;

            free(cutCoverLifted);
            free(solutionTemp);
        }
        // printf("antes do create: %d\n", binaryRows->numberConstraints);
        binaryRows = createCutsCover(cutsCoverSolution, binaryRows, copyRowsBinary, i, c_AuxSolution, precision);

        // printf("depois do create: %d\n", binaryRows->numberConstraints);
        freeStrConstraintsReal(cutsCoverSolution);
        free(poolSolution);
        currentTime = omp_get_wtime() - startT;
        if (currentTime > timeLeft)
        {
            break;
        }
    }
    //free(intOrFloat);
    //free(binaryRows);
    freeStrConstraintsReal(copyRowsBinary);
    // printf("depois de Tudo 2: %d\n", binaryRows->numberConstraints);
    return binaryRows;
    //retornar a solução
}

constraintsReal *returnVariablesOriginals(constraintsReal *constraintsOriginal, int *convertVector, int precision, int nVariablesInitial)
{
    //printf("%d - nVariablesInitial\n", nVariablesInitial);
    constraintsReal *newConstraints = AllocStrConstraintsReal(constraintsOriginal->cont, constraintsOriginal->numberConstraints, nVariablesInitial);
    int i, j, el;
    double rhs;
    for (i = 0; i < nVariablesInitial; i++)
    {
        newConstraints->xAsterisc[i] = constraintsOriginal->xAsterisc[i];
        newConstraints->ub[i] = constraintsOriginal->ub[i];
        newConstraints->lb[i] = constraintsOriginal->lb[i];
        newConstraints->typeVariables[i] = constraintsOriginal->typeVariables[i];
    }
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        rhs = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            el = constraintsOriginal->Elements[j];
            //printf("valor de el: %d\n",el);
            if (el >= nVariablesInitial)
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j] * (-1);
                rhs -= constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = convertVector[el - nVariablesInitial];
                newConstraints->xAsterisc[newConstraints->Elements[j]] = 1 - constraintsOriginal->xAsterisc[el]; //erro aqui
            }
            else
            {
                newConstraints->Coefficients[j] = constraintsOriginal->Coefficients[j];
                newConstraints->Elements[j] = constraintsOriginal->Elements[j];
            }
        }

        newConstraints->rightSide[i] = rhs;
        newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    }
    newConstraints->ElementsConstraints[i] = constraintsOriginal->ElementsConstraints[i];
    freeStrConstraintsReal(constraintsOriginal);
    return newConstraints;
}

constraintsReal *convertBinaryOfOriginalConstraints(constraintsReal *constraintsOriginal, constraintsReal *constraintsBinary, int nInitialBinary)
{
    int sizeNew = constraintsOriginal->numberConstraints + (constraintsBinary->numberConstraints - nInitialBinary);
    int contNew = constraintsOriginal->cont;
    int i, j, jTemp;
    for (i = nInitialBinary; i < constraintsBinary->numberConstraints; i++)
    {
        contNew += constraintsBinary->ElementsConstraints[i + 1] - constraintsBinary->ElementsConstraints[i];
    }
    if (constraintsBinary->numberConstraints - nInitialBinary == 0)
    {
        return constraintsOriginal;
    }
    // printf("Mais um teste: %d %d\n", constraintsBinary->numberConstraints, nInitialBinary);
    constraintsReal *newConstraintsOriginal = AllocStrConstraintsReal(contNew, sizeNew, constraintsOriginal->numberVariables);
    // printf("szConstraints %d - %d contNew %d sizeNew %d\n nInitial: %d",constraintsOriginal->numberConstraints,newConstraintsOriginal->numberConstraints, contNew,sizeNew, nInitialBinary);
    for (j = 0; j < constraintsOriginal->numberVariables; j++)
    {
        newConstraintsOriginal->xAsterisc[j] = constraintsOriginal->xAsterisc[j];
        newConstraintsOriginal->ub[j] = constraintsOriginal->ub[j];
        newConstraintsOriginal->lb[j] = constraintsOriginal->lb[j];
        newConstraintsOriginal->typeVariables[j] = constraintsOriginal->typeVariables[j];
    }
    newConstraintsOriginal->ElementsConstraints[0] = constraintsOriginal->ElementsConstraints[0];
    for (i = 0; i < constraintsOriginal->numberConstraints; i++)
    {
        newConstraintsOriginal->ElementsConstraints[i + 1] = constraintsOriginal->ElementsConstraints[i + 1];
        newConstraintsOriginal->rightSide[i] = constraintsOriginal->rightSide[i];
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            newConstraintsOriginal->Coefficients[j] = constraintsOriginal->Coefficients[j];
            newConstraintsOriginal->Elements[j] = constraintsOriginal->Elements[j];
        }
    }
    contNew = newConstraintsOriginal->ElementsConstraints[i];
    jTemp = constraintsOriginal->numberConstraints;
    for (i = nInitialBinary; i < constraintsBinary->numberConstraints; i++)
    {
        for (j = constraintsBinary->ElementsConstraints[i]; j < constraintsBinary->ElementsConstraints[i + 1]; j++)
        {
            newConstraintsOriginal->Coefficients[contNew] = constraintsBinary->Coefficients[j];
            newConstraintsOriginal->Elements[contNew] = constraintsBinary->Elements[j];
            contNew++;
        }
        newConstraintsOriginal->ElementsConstraints[jTemp + 1] = contNew;
        newConstraintsOriginal->rightSide[jTemp] = constraintsBinary->rightSide[i];
        jTemp++;
    }

    freeStrConstraintsReal(constraintsOriginal);
    return newConstraintsOriginal;
}

std::vector<std::string> renamedNameConstraints(std::vector<std::string> nameConstraints, int typeContraints, int szNewConstraints, int szCuts, int lastCut)
{
    int i;

    std::vector<std::string> newNameConstraints;
    for (i = 0; i < szNewConstraints - szCuts; i++)
    {
        newNameConstraints.push_back(nameConstraints[i]);
    }
    for (i = szNewConstraints - szCuts; i < szNewConstraints; i++)
    {
        if (typeContraints == 1)
        {
            std::string name;
            name = "CG1(";
            name.append(std::to_string(lastCut));
            name.append(")");
            //std::cout<<name<<std::endl;
            lastCut++;
            newNameConstraints.push_back(name);
        }
        if (typeContraints == 2)
        {
            std::string name;
            name = "CG2(";
            name.append(std::to_string(lastCut));
            name.append(")");
            lastCut++;
            newNameConstraints.push_back(name);
        }
        if (typeContraints == 3)
        {
            std::string name;
            name = "CK(";
            name.append(std::to_string(lastCut));
            name.append(")");
            lastCut++;
            //std::cout<<name<<std::endl;
            newNameConstraints.push_back(name);
            //std::cout<<name<<newNameConstraints[i]<<std::endl;
        }
    }
    //getchar();
    return newNameConstraints;
}

//constraintsReal *LiftWithSA(constraintsReal *constraintsInitial /*,std::vector<std::string> nameConstraints, std::vector<std::string> nameVariables*/)
constraintsReal *LiftWithSA(constraintsReal *constraintsInitial, 
                                double tempInitial, double fatorResf, double pv1,double pv2, double pv3, double pv4, double pv5,
                                double timeConstraints, int fm, double timeLeft, int typeLift, int typeSolutionInitial,
                                int precision, int minimal, int szPoolCutsMaxCC)
{

    // float tempInitial = 1.84371407, fatorResf = 0.99996506;                                         //Parameters of SA
    // float pv1 = 0.70270979, pv2 = 0.52492615, pv3 = 0.30411066, pv4 = 0.11879197, pv5 = 0.91005716; // Probabiity of neighboord -LACH-SA-SCHC
    // float sPr = pv1 + pv2 + pv3 + pv4 + pv5;
    // if (sPr != 0.0)
    // {
    //     pv1 = pv1 / sPr;
    //     pv2 = pv2 / sPr;
    //     pv3 = pv3 / sPr;
    //     pv4 = pv4 / sPr;
    //     pv5 = pv5 / sPr;
    // }
    // else
    // {
    //     pv1 = 0.2;
    //     pv2 = 0.2;
    //     pv3 = 0.2;
    //     pv4 = 0.2;
    //     pv5 = 0.2;
    // }
    // int szPoolCutsMaxCC = 200000; //max number of cuts
    // int precision = 1000;
    // int minimal = 0;             // 0 - no minimal 1- minimal
    // int typeSolutionInitial = 0; // type of prioritary greedy
    // int typeLift = 0;            //  0 - Adam Lethford 1 - Bala1s
    // double timeConstraints = 0.5;
    // int fm = 355;
    // //double timeLeft = timeConstraints * constraintsInitial->numberConstraints;
    // double timeLeft = 60;
    int nRowsOg = constraintsInitial->numberConstraints;
    constraintsReal *binaryConstraints = prepareConstraintsOfLifting(constraintsInitial, precision);
    // for(int i=0;i<binaryConstraints->numberConstraints;i++){
    //   for (int j = binaryConstraints->ElementsConstraints[i];j<binaryConstraints->ElementsConstraints[i+1];j++){
    //     printf("%lf x %d = %lf + ", binaryConstraints->Coefficients[j], binaryConstraints->Elements[j], binaryConstraints->xAsterisc[ binaryConstraints->Elements[j]]);
    //   }
    //   printf("<= %lf\n", binaryConstraints->rightSide[i]);
    // }
    int nRowsInitial = binaryConstraints->numberConstraints;
    int nColsInitial = binaryConstraints->numberVariables;
    int *convertVariables = (int *)malloc(sizeof(int) * binaryConstraints->cont);
    binaryConstraints = removeNegativeCoefficientsAndSort(binaryConstraints, convertVariables);
    int numberAuxRows = binaryConstraints->numberConstraints;
    //printf("number aux rows!!:%d %d %d\n", numberAuxRows, nRowsInitial, nRowsOg);
    double violationFinal = 0;
    binaryConstraints = runCCwithSA(binaryConstraints, precision, szPoolCutsMaxCC, tempInitial, fatorResf, typeLift, minimal, pv1, pv2, pv3, pv4, &violationFinal, typeSolutionInitial, timeLeft, fm, timeConstraints);
    numberAuxRows = binaryConstraints->numberConstraints - numberAuxRows;
    // printf("aqui vai mais um: %d %d\n", nRowsInitial,numberAuxRows);
    // printf("teste: %d %d %d\n",numberAuxRows, binaryConstraints->numberConstraints,binaryConstraints->numberVariables);
    binaryConstraints = returnVariablesOriginals(binaryConstraints, convertVariables, precision, nColsInitial);
    constraintsInitial = convertBinaryOfOriginalConstraints(constraintsInitial, binaryConstraints, nRowsInitial);
    // if (numberAuxRows > 0)
    // {
    //   nameConstraints = renamedNameConstraints(nameConstraints, 3, constraintsInitial->numberConstraints, numberAuxRows, nRowsOg);
    // }
    free(convertVariables);
    freeStrConstraintsReal(binaryConstraints);
    return constraintsInitial;
}

constraintsReal *useMatrixConflict(constraintsReal *constraintsInitial, std::vector<std::vector<int>> matrixConflict, std::map<int, int> &varXclique)
{
    int nVariablesOriginal = constraintsInitial->numberVariables;
    int szCliques = matrixConflict.size();
    int nVariablesNew = nVariablesOriginal + szCliques;
    int constraint;
    int i, j, k;
    std::vector<int> szConstraints;
    int *valitedConstraints = (int *)malloc(sizeof(int) * constraintsInitial->numberConstraints);
    for (i = 0; i < constraintsInitial->numberConstraints; i++)
    {
        valitedConstraints[i] = 1;
        int szAux = constraintsInitial->ElementsConstraints[i + 1] - constraintsInitial->ElementsConstraints[i];
        szConstraints.push_back(szAux);
    }
    int nContAux = 0, nConsAux = 0;
    for (i = 0; i < szCliques; i++)
    {
        double xAst = 0;
        int aux = matrixConflict[i].size();
        for (j = 0; j < aux - 1; j++)
        {
            int el = matrixConflict[i][j];
            //printf("%f x%d = %f +", constraintsInitial->Coefficients[el], constraintsInitial->Elements[el], constraintsInitial->xAsterisc[constraintsInitial->Elements[el]]);
        }
        constraint = matrixConflict[i][aux - 1];
        szConstraints[constraint] -= (aux - 2);
        valitedConstraints[constraint] = 2;
        //printf("\nconstraints: %d\n",constraint);
    }
    for (i = 0; i < constraintsInitial->numberConstraints; i++)
    {
        if (szConstraints[i] > 1)
        {
            nContAux += szConstraints[i];
            nConsAux++;
        }
        else
        {
            valitedConstraints[i] = 0;
        }
    }
    constraintsReal *newConstraints = AllocStrConstraintsReal(nContAux, nConsAux, nVariablesNew);
    newConstraints->ElementsConstraints[0] = 0;
    for (i = 0; i < constraintsInitial->numberVariables; i++)
    {
        newConstraints->xAsterisc[i] = constraintsInitial->xAsterisc[i];
        newConstraints->ub[i] = constraintsInitial->ub[i];
        newConstraints->lb[i] = constraintsInitial->lb[i];
        newConstraints->typeVariables[i] = constraintsInitial->typeVariables[i];
    }
    int posAux = 0, constAux = 0;
    for (i = 0; i < constraintsInitial->numberConstraints; i++)
    {
        if (valitedConstraints[i] == 1)
        {
            for (k = constraintsInitial->ElementsConstraints[i]; k < constraintsInitial->ElementsConstraints[i + 1]; k++)
            {
                newConstraints->Coefficients[posAux] = constraintsInitial->Coefficients[k];
                newConstraints->Elements[posAux] = constraintsInitial->Elements[k];
                posAux++;
            }
            newConstraints->ElementsConstraints[constAux + 1] = posAux;
            newConstraints->rightSide[constAux] = constraintsInitial->rightSide[i];
            constAux++;
        }
        else if (valitedConstraints[i] == 2)
        {
            int flag = 0;
            std::set<int> position;
            std::vector<std::pair<double, int>> inserirNova;
            for (k = constraintsInitial->ElementsConstraints[i]; k < constraintsInitial->ElementsConstraints[i + 1]; k++)
            {
                position.insert(k);
            }
            for (j = 0; j < szCliques; j++)
            {
                int aux = matrixConflict[j].size();
                if (i == matrixConflict[j][aux - 1])
                {

                    double xAstAux = 0;
                    int cliqueAux = j; //constraintsInitial->Coefficients[ matrixConflict[j][0]];
                    for (k = 0; k < aux - 1; k++)
                    {
                        xAstAux += constraintsInitial->xAsterisc[constraintsInitial->Elements[matrixConflict[j][k]]];
                        std::set<int>::iterator it;
                        it = position.find(matrixConflict[j][k]);
                        position.erase(it);
                    }
                    flag = 1;
                    std::pair<double, int> aux(xAstAux, cliqueAux);
                    inserirNova.push_back(aux);
                }
                else if (flag == 1)
                {
                    break;
                }
            }
            for (std::set<int>::iterator it = position.begin(); it != position.end(); ++it)
            {
                newConstraints->Coefficients[posAux] = constraintsInitial->Coefficients[*it];
                newConstraints->Elements[posAux] = constraintsInitial->Elements[*it];
                posAux++;
            }
            for (k = 0; k < inserirNova.size(); k++)
            {
                newConstraints->Coefficients[posAux] = constraintsInitial->Coefficients[matrixConflict[inserirNova[k].second][0]];
                newConstraints->Elements[posAux] = nVariablesOriginal;
                varXclique.insert(std::pair<int, int>(nVariablesOriginal, inserirNova[k].second));
                newConstraints->xAsterisc[nVariablesOriginal] = inserirNova[k].first;
                newConstraints->ub[nVariablesOriginal] = 1;
                newConstraints->lb[nVariablesOriginal] = 0;
                newConstraints->typeVariables[nVariablesOriginal] = 0;
                nVariablesOriginal++;
                posAux++;
            }

            newConstraints->ElementsConstraints[constAux + 1] = posAux;
            newConstraints->rightSide[constAux] = constraintsInitial->rightSide[i];
            constAux++;
        }
    }
    free(valitedConstraints);
    return newConstraints;
}

constraintsReal *returnOriginalConflict(constraintsReal *cOriginal, constraintsReal *cNew, std::vector<std::vector<int>> matrixConflict, std::map<int, int> varXclique, int nCuts)
{
    int nCOriginal = cOriginal->numberConstraints;
    int nVOriginal = cOriginal->numberVariables;
    int i, j, k, cont = 0;
    for (i = cNew->numberConstraints - nCuts; i < cNew->numberConstraints; i++)
    {
        for (j = cNew->ElementsConstraints[i]; j < cNew->ElementsConstraints[i + 1]; j++)
        {
            if (cNew->Elements[j] >= nVOriginal)
            {

                int clique = varXclique[cNew->Elements[j]];
                cont += matrixConflict[clique].size();
            }
            else
            {
                cont++;
            }
        }
    }
    constraintsReal *outputConstraints = AllocStrConstraintsReal(cOriginal->cont + cont, nCOriginal + nCuts, nVOriginal);
    outputConstraints->ElementsConstraints[0] = 0;
    for (i = 0; i < nVOriginal; i++)
    {
        outputConstraints->xAsterisc[i] = cOriginal->xAsterisc[i];
        outputConstraints->ub[i] = cOriginal->ub[i];
        outputConstraints->lb[i] = cOriginal->lb[i];
        outputConstraints->typeVariables[i] = cOriginal->typeVariables[i];
    }
    for (i = 0; i < nCOriginal; i++)
    {
        for (j = cOriginal->ElementsConstraints[i]; j < cOriginal->ElementsConstraints[i + 1]; j++)
        {
            outputConstraints->Coefficients[j] = cOriginal->Coefficients[j];
            outputConstraints->Elements[j] = cOriginal->Elements[j];
        }
        outputConstraints->ElementsConstraints[i + 1] = cOriginal->ElementsConstraints[i + 1];
        outputConstraints->rightSide[i] = cOriginal->rightSide[i];
    }
    cont = cOriginal->cont;
    int cAux = cOriginal->numberConstraints;
    //printf("NORRRR: %d\n", nVOriginal);
    for (i = cNew->numberConstraints - nCuts; i < cNew->numberConstraints; i++)
    {
        for (j = cNew->ElementsConstraints[i]; j < cNew->ElementsConstraints[i + 1]; j++)
        {

            if (cNew->Elements[j] >= nVOriginal)
            {
                int clique = varXclique[cNew->Elements[j]];
                //printf("%d %d\n",clique,  cNew->Elements[j]);
                for (k = 0; k < matrixConflict[clique].size() - 1; k++)
                {
                    int posOriginal = matrixConflict[clique][k];
                    int el = cOriginal->Elements[posOriginal];
                    // printf("%d ", el);
                    outputConstraints->Coefficients[cont] = cNew->Coefficients[j];
                    outputConstraints->Elements[cont] = cOriginal->Elements[posOriginal];
                    cont++;
                }
            }
            else
            {
                outputConstraints->Coefficients[cont] = cNew->Coefficients[j];
                outputConstraints->Elements[cont] = cNew->Elements[j];
                //printf("%d ", cNew->Elements[j]);
                cont++;
            }
        }
        outputConstraints->ElementsConstraints[cAux + 1] = cont;
        outputConstraints->rightSide[cAux] = cNew->rightSide[i];
        cAux++;
    }
    freeStrConstraintsReal(cOriginal);
    freeStrConstraintsReal(cNew);
    return outputConstraints;
}

constraintsReal *readPseudoInstance(const char *nameFile, std::map<int, std::set<int>> &conflict)
{

    int i, cont;
    double coef, xAst, rhs, lb, ub;
    int el, cf1, cf2, tp;
    std::map<int, int> tempEl;
    // printf("%s\n", nameFile);
    // getchar();
    FILE *fin = fopen(nameFile, "r");
    fscanf(fin, "%d\n", &cont);
    constraintsReal *cOriginal = AllocStrConstraintsReal(cont, 1, cont);
    cOriginal->ElementsConstraints[0] = 0;
    for (i = 0; i < cont; i++)
    {
        fscanf(fin, "%lf %d %lf %d %lf %lf \n", &coef, &el, &xAst, &tp, &lb, &ub);
        cOriginal->Coefficients[i] = coef;
        cOriginal->Elements[i] = i;
        cOriginal->xAsterisc[i] = xAst;
        cOriginal->ub[i] = ub;
        cOriginal->lb[i] = lb;
        cOriginal->typeVariables[i] = tp;
        tempEl[el] = i;
    }
    cOriginal->ElementsConstraints[1] = i;
    fscanf(fin, "%lf\n", &rhs);
    cOriginal->rightSide[0] = rhs;
    while (!feof(fin))
    {
        fscanf(fin, "%d %d\n", &cf1, &cf2);
        conflict[tempEl[cf1]].insert(tempEl[cf2]);
        conflict[tempEl[cf2]].insert(tempEl[cf1]);
    }
    fclose(fin);
    return cOriginal;
}

bool conflicting(std::map<int, std::set<int>> conflictMapSet, int el1, int el2)
{
    if (conflictMapSet[el1].size() == 0)
    {
        return 0;
    }
    else
    {
        for (std::map<int, std::set<int>>::iterator it = conflictMapSet.begin(); it != conflictMapSet.end(); ++it)
        {
            if (el1 == it->first)
            {
                for (std::set<int>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2)
                {
                   if(*it2==el2){
                       return 1;
                   }
                }
                return 0;
            }
        }
    }
}

int *vectorNonRepeteadNonDominated(constraintsReal *constraintsOriginal, int nConstraintsInitial)
{
    int i, j, k = 0, qntRepeat = 0, qntDominated = 0;
   // printf("ok guy: %d\n", nConstraintsInitial);
    int *isRepeat = (int *)malloc(sizeof(int) * constraintsOriginal->numberConstraints);
    for (i = 0; i < nConstraintsInitial; i++)
    {
        isRepeat[i] = 0;
    }
    //int *v_aux = (int *)malloc(sizeof(int) * (constraintsOriginal->numberVariables+1) );
    std::vector <double> v_aux;
    for (i = nConstraintsInitial; i < constraintsOriginal->numberConstraints; i++)
    {   
        k = 0;
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            //v_aux[k] = constraintsOriginal->Coefficients[j];
            v_aux.push_back(constraintsOriginal->Coefficients[j]);
            k++;
        }
        //v_aux[k] = constraintsOriginal->rightSide[i];
        v_aux.push_back(constraintsOriginal->rightSide[i]);
        int mdc = cutMaxDivisorCommonVector(v_aux, k);
        for (j = constraintsOriginal->ElementsConstraints[i]; j < constraintsOriginal->ElementsConstraints[i + 1]; j++)
        {
            constraintsOriginal->Coefficients[j] = constraintsOriginal->Coefficients[j] / mdc;
            // printf("%lf x%d +",constraintsOriginal->Coefficients[j], constraintsOriginal->Elements[j] );
        }
        constraintsOriginal->rightSide[i] = constraintsOriginal->rightSide[i] / mdc;
        // printf("<= %lf\n",constraintsOriginal->rightSide[i]);
        for (j = 0; j < i; j++)
        {
            isRepeat[i] = verifyRepeatCuts(constraintsOriginal, j, i);
            qntRepeat += isRepeat[i];
            if (isRepeat[i] == 1)
            {
                break;
            }
        }
    }
    for (i = nConstraintsInitial; i < constraintsOriginal->numberConstraints; i++)
    {  
        if(isRepeat[i]==0){
            for (j = 0; j < constraintsOriginal->numberConstraints; j++)
            {  
                if(i!=j){
                    isRepeat[i] = verifyDominanceCG(constraintsOriginal, j, i);
                    if(isRepeat[i]==1){
                        qntDominated++;
                        break;
                    }
                }
            }
        }
    }


    // printf("ncRpt: %d\t", qntRepeat);
    // printf("ncDominated: %d\t", qntDominated);
    // //free(v_aux);
    return isRepeat;
}
    
int verifyRepeatCuts(constraintsReal *constraintsOriginal, int cutOriginal, int cutCreate)
{
    int i, j, aux = 1;
    int szOri = constraintsOriginal->ElementsConstraints[cutOriginal + 1] - constraintsOriginal->ElementsConstraints[cutOriginal];
    int szCre = constraintsOriginal->ElementsConstraints[cutCreate + 1] - constraintsOriginal->ElementsConstraints[cutCreate];
    if ((szOri != szCre) || (constraintsOriginal->rightSide[cutOriginal] != constraintsOriginal->rightSide[cutCreate]))
    {
        return 0;
    }
    for (i = constraintsOriginal->ElementsConstraints[cutOriginal]; i < constraintsOriginal->ElementsConstraints[cutOriginal + 1]; i++)
    {
        aux = 0;
        for (j = constraintsOriginal->ElementsConstraints[cutCreate]; j < constraintsOriginal->ElementsConstraints[cutCreate + 1]; j++)
        {
            if ((constraintsOriginal->Coefficients[i] == constraintsOriginal->Coefficients[j]) && (constraintsOriginal->Elements[i] == constraintsOriginal->Elements[j]))
            {
                aux = 1;
            }
        }
        if (aux == 0)
        {
            return 0;
        }
    }
    return 1;
}

    

/*verify dominance between two vector*/
/*c1 dominate c2*/
int verifyDominanceCG(constraintsReal *constraintsInitial, int c1, int c2)
{
    int *v1 = (int*)malloc(sizeof(int)*constraintsInitial->numberVariables);
    int *v2 = (int*)malloc(sizeof(int)*constraintsInitial->numberVariables);
    int i, cont = 0, el;
    for(i= 0; i<constraintsInitial->numberVariables;i++){
        v1[i] = 0;
        v2[i] = 0;
    }
    for (i = constraintsInitial->ElementsConstraints[c1]; i<constraintsInitial->ElementsConstraints[c1+1];i++){
        el = constraintsInitial->Elements[i];
        v1[el] = constraintsInitial->Coefficients[i];
    }
    for (i = constraintsInitial->ElementsConstraints[c2]; i<constraintsInitial->ElementsConstraints[c2+1];i++){
        el = constraintsInitial->Elements[i];
        v2[el] = constraintsInitial->Coefficients[i];
    }
    int rhs1 = constraintsInitial->rightSide[c1];
    int rhs2 = constraintsInitial->rightSide[c2];
    int sz = constraintsInitial->numberVariables;
    

    for(i=0; i<sz; i++)
    {
        v1[i] = v1[i]*rhs2;
        v2[i] = v2[i]*rhs1;
    }
    rhs2 = rhs2*rhs1;
    rhs1 = rhs2;
    //int dominance = 1;
    for(i=0; i<sz; i++)
    {
        if(v1[i]<v2[i])
        {   
            free(v1);
            free(v2);
            return 0;
        }
        else if(v1[i]>v2[i])
        {
            cont++;
        }
    }
    free(v1);
    free(v2);
    if(cont>0)
        return 1;
    else
        return 0;
}



double countViolation(constraintsReal *constraintsInitial, int *constraintsValited, int nConstraints){
    double violation = 0, vAux = 0;
    int i, j, el;
    for(i=nConstraints;i<constraintsInitial->numberVariables;i++){
        if(constraintsValited[i]==0){
            vAux = 0;
            for(j = constraintsInitial->ElementsConstraints[i];j<constraintsInitial->ElementsConstraints[i+1];j++){
                el  = constraintsInitial->Elements[j];
                vAux += constraintsInitial->Coefficients[j]*constraintsInitial->xAsterisc[el];
                //printf("%lf %lf \n",constraintsInitial->Coefficients[j], constraintsInitial->xAsterisc[el]);
                
            }
            if(vAux - constraintsInitial->rightSide[i]>0 ){
                //printf("vAux: %lf rhs: %lf nVariables: %d\n",vAux, constraintsInitial->rightSide[i], constraintsInitial->numberVariables);
                violation += vAux - constraintsInitial->rightSide[i];
            }
        }
    }
    return violation;
}