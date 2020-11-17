#include "functions.hpp"

int main(int argc, const char *argv[])
{
    char nameFileInstance[255] = "";
    int i = 0;
    double tempInitial = 1.84371407, fatorResf = 0.99996506;  
    double pv1 = 0.70270979, pv2 = 0.52492615, pv3 = 0.30411066, pv4 = 0.11879197, pv5 = 0.91005716;
    int szPoolCutsMaxCC = 200000; //max number of cuts
    int precision = 1000;
    int minimal = 0;             // 0 - no minimal 1- minimal
    int typeSolutionInitial = 0; // type of prioritary greedy
    int typeLift = 0;            //  0 - Adam Lethford 1 - Bala1s
    double timeConstraints = 0.5;
    int fm = 355;
    //double timeLeft = timeConstraints * constraintsInitial->numberConstraints;
    double timeLeft = 60;
    int nIterationConflict = 100;
    double weightTime = 1;
    for (i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-f") == 0)
        { //Final source
            strcat(nameFileInstance, argv[i + 1]);
        }
        // if (strcmp(argv[i], "-t") == 0)
        // {
        //     // time Limit
        //     timeMax = atof(argv[i + 1]);
        // }
        if (strcmp(argv[i], "-tc") == 0)
        {
            // time Limit per Constraints
            timeConstraints = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-p") == 0)
        {
            //precision of integer number
            precision = atoi(argv[i + 1]);
        }
        if (strcmp(argv[i], "-l") == 0)
        { //type lifting 0 -  Adam e 1 - Ballas
            typeLift = atoi(argv[i + 1]);
        }
        if (strcmp(argv[i], "-m") == 0)
        { //minimal cover set- 0 no minimal and 1 minimal
            minimal = atoi(argv[i + 1]);
        }
        if (strcmp(argv[i], "-sp") == 0)
        { //sizePool
            szPoolCutsMaxCC = atoi(argv[i + 1]);
        }
        // if (strcmp(argv[i], "-i") == 0)
        // { //number iteration diverses heuristics
        //     nIterationCC = atoi(argv[i + 1]);
        // }
        // if (strcmp(argv[i], "-a") == 0)
        // { //alpha GRASP
        //     alpha = atof(argv[i + 1]);
        // }
        if (strcmp(argv[i], "-ti") == 0)
        { //initial temperature SA
            tempInitial = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-fr") == 0)
        { //factor reduce SA
            fatorResf = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-pv1") == 0)
        { //probality neighborhood 5
            pv1 = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-pv2") == 0)
        { //probality neighborhood 2
            pv2 = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-pv3") == 0)
        { //probality neighborhood 3
            pv3 = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-pv4") == 0)
        { //probality neighborhood 4
            pv4 = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-pv5") == 0)
        { //probality neighborhood 5
            pv5 = atof(argv[i + 1]);
        }
        // if (strcmp(argv[i], "-LFA") == 0)
        // { //LFA LAHC
        //     LFA = atoi(argv[i + 1]);
        // }
        if (strcmp(argv[i], "-tsi") == 0)
        { //Type initial Solution
            typeSolutionInitial = atoi(argv[i + 1]);
        }
        // if (strcmp(argv[i], "-gap") == 0)
        // { //Gap Verify
        //     gap = atoi(argv[i + 1]);
        // }
        // if (strcmp(argv[i], "-sc") == 0)
        // { //set Cut
        //     setCut = atoi(argv[i + 1]);
        // }
        // if (strcmp(argv[i], "-lc") == 0)
        // { //List Counting SCHC
        //     lc = atoi(argv[i + 1]);
        // }
        if (strcmp(argv[i], "-fm") == 0)
        { //factor mutiplier no improvement
            fm = atoi(argv[i + 1]);
        }
        if (strcmp(argv[i], "-tl") == 0)
        { //time max SA
            timeLeft = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "-nCf") == 0)
        { //time max SA
            nIterationConflict = atoi(argv[i + 1]);
        }
        if (strcmp(argv[i], "-wt") == 0)
        { //time max SA
            weightTime = atof(argv[i + 1]);
        }
    }
    std::map<int, std::set<int>> conflictMapSet;
    // printf("name %s\n", nameFileInstance);
    constraintsReal *constraintsInitial = readPseudoInstance(nameFileInstance, conflictMapSet);
    //  for (i = 0; i < constraintsInitial->cont; i++)
    //  {
    //      int el = constraintsInitial->Elements[i];
    //      printf("%lf x%d = %lf + ", constraintsInitial->Coefficients[i], el, constraintsInitial->xAsterisc[el]);
    //  }
    //  printf("<= %lf\n", constraintsInitial->rightSide[0]);

    //int nIterationConflict = 100;
    int nVariablesConflict = constraintsInitial->numberVariables;
    std::vector<std::vector<int>> matrixConflict;
    for (int i = 0; i < constraintsInitial->numberConstraints; i++)
    {
        std::map<float, std::set<int>> mapSetCoef;
        std::map<float, std::set<int>> mapSetPosition;
        std::map<int, std::set<int>> mapSetConstraintsVariables;
        for (int j = constraintsInitial->ElementsConstraints[i]; j < constraintsInitial->ElementsConstraints[i + 1]; j++)
        {
            float coeftemp = constraintsInitial->Coefficients[j];
            int el = constraintsInitial->Elements[j];
            (mapSetCoef[coeftemp]).insert(el);
            (mapSetPosition[coeftemp]).insert(j);
        }
        for (std::map<float, std::set<int>>::iterator it = mapSetCoef.begin(); it != mapSetCoef.end(); ++it)
        {
            //std::cout<< it->first<<"("<<(it->second).size()<<"):";
            std::vector<int> auxSet;
            std::vector<int> auxSetPosition;
            std::vector<int> auxSetConflict;
            for (std::set<int>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2)
            {

                auxSet.push_back(*it2);

                //std::cout << ' ' << *it2;
            }
            for (std::set<int>::iterator it2 = (mapSetPosition[it->first]).begin(); it2 != (mapSetPosition[it->first]).end(); ++it2)
            {

                auxSetPosition.push_back(*it2);
                //std::cout << ' ' << *it2;
            }

            int iteConf = 0;
            int flag = 0;
            while ((iteConf < nIterationConflict) && (auxSet.size() > 1))
            {
                if (flag == 1)
                {
                    iteConf++;
                }
                while ((auxSet.size() > 1) && (flag == 0) && (iteConf < nIterationConflict))
                {
                    int v1, v2;
                    v1 = rand() % auxSet.size();
                    do
                    {
                        v2 = rand() % auxSet.size();
                    } while (v1 == v2);
                    if (conflicting(conflictMapSet, auxSet[v1], auxSet[v2]))
                    {

                        auxSetConflict.push_back(auxSetPosition[v1]);
                        auxSetConflict.push_back(auxSetPosition[v2]);
                        if (v1 < v2)
                        {
                            auxSet.erase(auxSet.begin() + v1);
                            auxSet.erase(auxSet.begin() + (v2 - 1));
                            auxSetPosition.erase(auxSetPosition.begin() + v1);
                            auxSetPosition.erase(auxSetPosition.begin() + (v2 - 1));
                        }
                        else
                        {
                            auxSet.erase(auxSet.begin() + v1);
                            auxSet.erase(auxSet.begin() + v2);
                            auxSetPosition.erase(auxSetPosition.begin() + v1);
                            auxSetPosition.erase(auxSetPosition.begin() + v2);
                        }

                        flag = 1;
                    }
                    iteConf++;
                }
                if ((auxSet.size() < 2) || (iteConf >= nIterationConflict))
                {
                    break;
                }
                int v = rand() % auxSet.size();
                bool checkConf = true;
                for (int j = 0; j < auxSetConflict.size(); j++)
                {
                    int el = constraintsInitial->Elements[auxSetConflict[j]];
                    if (conflicting(conflictMapSet, auxSet[v], el) == 0)
                    {
                        checkConf = false;
                        break;
                    }
                }
                if (checkConf = true)
                {
                    auxSetConflict.push_back(auxSetPosition[v]);
                    auxSet.erase(auxSet.begin() + v);
                    auxSetPosition.erase(auxSetPosition.begin() + v);
                }
            }
            if (auxSetConflict.size() >= 2)
            {
                auxSetConflict.push_back(i); //última posição sera a restrição
                matrixConflict.push_back(auxSetConflict);
            }
        }
    }

    int nCuts;
    double sPr = pv1 + pv2 + pv3 + pv4 + pv5;
    if (sPr != 0.0)
    {
        pv1 = pv1 / sPr;
        pv2 = pv2 / sPr;
        pv3 = pv3 / sPr;
        pv4 = pv4 / sPr;
        pv5 = pv5 / sPr;
    }
    else
    {
        pv1 = 0.2;
        pv2 = 0.2;
        pv3 = 0.2;
        pv4 = 0.2;
        pv5 = 0.2;
    }



    double startT = omp_get_wtime();
    
    if (!matrixConflict.empty())
    {
        std::map<int, int> varXclique;
        constraintsReal *newConstraints = useMatrixConflict(constraintsInitial, matrixConflict, varXclique);
        nCuts = newConstraints->numberConstraints;
        newConstraints = LiftWithSA(newConstraints , 
                                tempInitial,fatorResf,pv1,pv2,pv3,pv4,pv5,
                                timeConstraints, fm,timeLeft,typeLift,typeSolutionInitial,
                                precision, minimal, szPoolCutsMaxCC);
        nCuts = newConstraints->numberConstraints - nCuts;
        constraintsInitial = returnOriginalConflict(constraintsInitial, newConstraints, matrixConflict, varXclique, nCuts);
    }
    else
    {
        nCuts = constraintsInitial->numberConstraints;
        constraintsInitial = LiftWithSA(constraintsInitial , 
                                tempInitial,fatorResf,pv1,pv2,pv3,pv4,pv5,
                                timeConstraints, fm,timeLeft,typeLift,typeSolutionInitial,
                                precision, minimal, szPoolCutsMaxCC);
        nCuts = constraintsInitial->numberConstraints - nCuts;
    }

    int *noRp = vectorNonRepeteadNonDominated(constraintsInitial, constraintsInitial->numberConstraints - nCuts);
    double currentTime = omp_get_wtime() - startT;
    double violationFinal = 0;
    if(nCuts>0){
     violationFinal = countViolation(constraintsInitial,noRp,constraintsInitial->numberConstraints - nCuts);
    }
    //printf("%s nCuts: %d time %lf\n",nameFileInstance, nCuts, currentTime);
    //printf("ViolationFinal: %lf\n",violationFinal);
    double fitness =(-1)*(violationFinal - (weightTime*currentTime));
    printf("Fiteness: %lf\n", fitness);
    free(noRp);
    freeStrConstraintsReal(constraintsInitial);
}
