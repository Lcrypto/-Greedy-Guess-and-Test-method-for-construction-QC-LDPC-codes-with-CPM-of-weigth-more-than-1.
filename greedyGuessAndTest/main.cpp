/*
Copyright(c) 2012, Ilya Vorobyev und Vasiliy Usatyuk
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include".\myLib\irregularLDPC.h"
#include<queue>


struct Tiii {
    int first, second, third;
    Tiii() {}
    Tiii(int _first, int _second, int _third) : first(_first), second(_second), third(_third) {}
    Tiii(const entry& e) : first(e.r), second(e.c), third(e.id) {}
    bool operator==(const Tiii& rhs) {
        return (first == rhs.first) && (second == rhs.second) && (third == rhs.third);
    }
};

struct Cycle {
    vector<Tiii> cycle;
    int free;
    Cycle(const vector<entry>& _cycle) {
        cycle.resize(_cycle.size());
        for (size_t i = 0; i < cycle.size(); ++i)
            cycle[i] = _cycle[i];
        free = 1;
        for (int i = 1; i < cycle.size(); ++i) {
            bool unique = 1;
            for (int j = 0; j < i; ++j) {
                if (cycle[i] == cycle[j]) {
                    unique = 0;
                }
            }
            if (unique)
                ++free;
        }
    }
};



bool pickNext(Tiii& next, const vector<vector<vector<int> > >& numberOfPossibleValues, const vector<vector<vector<int> > >& assigned, ll circulant) {
    ll minimum = circulant + 1;
    for (int r = 0; r < assigned.size(); ++r) {
        for (int c = 0; c < assigned[r].size(); ++c) {
            for (int id = 0; id < assigned[r][c].size(); ++id) {
                if (assigned[r][c][id])
                    continue;
                if (numberOfPossibleValues[r][c][id] < minimum) {
                    minimum = numberOfPossibleValues[r][c][id];
                    next = Tiii(r, c, id);
                }
            }
        }
    }
    return (minimum > 0);
}

ll getValue(ll value, const vector<int>& possibleValues) {
    int counter = -1;
    for (int i = 0; i < possibleValues.size(); ++i) {
        if (possibleValues[i])
            ++counter;
        if (counter == value)
            return i;
    }
}

void processCycle(const Cycle& cycle, const vector<vector<vector<int> > >& assigned, const vector<vector<vector<int> > >& a, 
                    vector<vector<vector<vector<int> > > >& possibleValues, vector<vector<vector<int> > >& numberOfPossibleValues, ll circulant) {
    Tiii last;
    for (int i = 0; i < cycle.cycle.size(); ++i) {
        Tiii cur = cycle.cycle[i];
        if (assigned[cur.first][cur.second][cur.third])
            continue;
        last = cur;
        break;
    }
    ll c = 0, b = 0;//cx=b
    for (int i = 0; i < cycle.cycle.size(); ++i) {
        if (last == cycle.cycle[i]) {
            if (i & 1)
                ++c;
            else
                --c;
        }
        else {
            if (i & 1)
                b -= a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third];
            else
                b += a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third];
        }
    }
    c = ((c % circulant) + circulant) % circulant;
    b = ((b % circulant) + circulant) % circulant;
    for (int i = 0; i < possibleValues[last.first][last.second][last.third].size(); ++i) {
        if (possibleValues[last.first][last.second][last.third][i] == 0)
            continue;
        if ((c * i) % circulant == b) {
            possibleValues[last.first][last.second][last.third][i] = 0;
            --numberOfPossibleValues[last.first][last.second][last.third];
        }
    }
}

void processCycles(queue<pii>& cyclesToProcess, const vector<vector<Cycle> >& cycles, const vector<vector<vector<int> > >& assigned, const vector<vector<vector<int> > >& a,
                    vector<vector<vector<vector<int> > > >& possibleValues, vector<vector<vector<int> > >& numberOfPossibleValues, ll circulant) {
    while (!cyclesToProcess.empty()) {
        pii cycleId = cyclesToProcess.front();
        cyclesToProcess.pop();
        processCycle(cycles[cycleId.first][cycleId.second], assigned, a, possibleValues, numberOfPossibleValues, circulant);
     }
}

bool gen(int checkNodes, int variableNodes, const vector<vector<int> >& protograph, ll circulant, ll targetGirth, vector<vector<vector<int> > >& a) {
    a.assign(checkNodes, vector<vector<int> >(variableNodes));
    vector<vector<vector<vector<int> > > > possibleValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
    vector<vector<vector<int> > > assigned(checkNodes, vector<vector<int> >(variableNodes));
    vector<vector<vector<int> > > numberOfPossibleValues(checkNodes, vector<vector<int> >(variableNodes));
    vector<vector<vector<vector<pii> > > > cyclesByEntry(checkNodes, vector<vector<vector<pii> > > (variableNodes));;

    ll numberOfValuesToAssign = 0;
    for (int r = 0; r < checkNodes; ++r) {
        for (int c = 0; c < variableNodes; ++c) {
            a[r][c].resize(protograph[r][c]);
            numberOfValuesToAssign += protograph[r][c];
            possibleValues[r][c].assign(protograph[r][c], vector<int>(circulant, 1));
            assigned[r][c].assign(protograph[r][c], 0);
            numberOfPossibleValues[r][c].assign(protograph[r][c], circulant);
            cyclesByEntry[r][c].resize(protograph[r][c]);
        }
    }
    vector<vector<Cycle> > cycles(targetGirth);
    queue<pii> cyclesToProcess;
    for (int girth = 4; girth < targetGirth; girth += 2) {
        CycleEnum enumerator(girth, protograph);
        if (!enumerator.init()) {
            continue;
        }
        do {
            cycles[girth].push_back(Cycle(enumerator.cycle));
            for (int i = 0; i < girth; ++i) {
                cyclesByEntry[enumerator.cycle[i].r][enumerator.cycle[i].c][enumerator.cycle[i].id].push_back(pii(girth, cycles[girth].size() - 1));
            }
            if (cycles[girth][cycles[girth].size() - 1].free == 1)
                cyclesToProcess.push(pii(girth, cycles[girth].size() - 1));
        } while (enumerator.next());
    }
    for (int r = 0; r < checkNodes; ++r) {
        for (int c = 0; c < variableNodes; ++c) {
            for (int id = 0; id < protograph[r][c]; ++id) {
                //cerr << r << "\t" << c << "\t" << id << endl;
                sort(cyclesByEntry[r][c][id].begin(), cyclesByEntry[r][c][id].end());
                cyclesByEntry[r][c][id].resize(unique(cyclesByEntry[r][c][id].begin(), cyclesByEntry[r][c][id].end()) - cyclesByEntry[r][c][id].begin());
            }
        }
    }
    Tiii next;
    for (int stepId = 0; stepId < numberOfValuesToAssign; ++stepId) {
        processCycles(cyclesToProcess, cycles, assigned, a, possibleValues, numberOfPossibleValues, circulant);
        if (!pickNext(next, numberOfPossibleValues, assigned, circulant))
            return 0;
        ll value = getRand(numberOfPossibleValues[next.first][next.second][next.third]);
        a[next.first][next.second][next.third] = getValue(value, possibleValues[next.first][next.second][next.third]);
        assigned[next.first][next.second][next.third] = 1;
        for (int i = 0; i < cyclesByEntry[next.first][next.second][next.third].size(); ++i) {
            pii cycleId = cyclesByEntry[next.first][next.second][next.third][i];
            --cycles[cycleId.first][cycleId.second].free;
            if (cycles[cycleId.first][cycleId.second].free == 1)
                cyclesToProcess.push(cycleId);
        }
    }
    return 1;
}





bool girthAtLeast6Enum(const vector<vector<vector<int> > >& a, const vector<vector<int> >& protograph, ll circulant) {
    CycleEnum cycle(4, protograph);
    if (!cycle.init())
        return 1;
    while (true) {
        int res = 0;
        for (int i = 0; i * 2 < 4; ++i) {
            res = res + circulant + a[cycle.cycle[2 * i].r][cycle.cycle[2 * i].c][cycle.cycle[2 * i].id] - a[cycle.cycle[2 * i + 1].r][cycle.cycle[2 * i + 1].c][cycle.cycle[2 * i + 1].id];
        }
        if (res % circulant == 0) {
            /*for (int i = 0; i < 4; ++i) {
            cout << cycle.cycle[i].r << " " << cycle.cycle[i].c << " " << cycle.cycle[i].id << endl;
            }*/
            return 0;

        }
        if (!cycle.next())
            return 1;
    }
}

void test() {
    time_t start = time(NULL);
    int variableNodes = 9;
    int checkNodes = 5;
    vector<vector<int> > protograph = {
        { 0, 0, 0, 1, 0, 0, 0, 1, 2 },
        { 0, 0, 0, 0, 0, 1, 1, 1, 2 },
        { 0, 1, 0, 0, 1, 1, 1, 0, 2 },
        { 0, 0, 1, 0, 1, 1, 0, 1, 2 },
        { 2, 1, 1, 1, 1, 0, 1, 0, 1 }};
    ll iterationCount = 0;
    ll circulant = 42;
    ll seed = 123;
    srand(seed);
    while (true) {
        ++iterationCount;
        vector<vector<vector<int> > > a;
        gen(checkNodes, variableNodes, protograph, circulant, 6, a);
        //cerr << getGirth(a, protograph, circulant) << endl;
        bool res1 = girthAtLeast6Manual(a, circulant);
        bool res2 = girthAtLeast6Enum(a, protograph, circulant);
        if (res1 != res2) {
            cerr << "ERROR" << endl;
            cout << res1 << endl;
            cout << res2 << endl;
            print(a);
            break;
        }
        if (iterationCount % 1000000 == 0) {
            cout << iterationCount / 1000000 << " million iterations\n";
            cout << time(NULL) - start << " seconds\n";
        }
    }

}





int main(int argc, char* argv[]) {
    bool validInput = 1;
    bool regular = 0;
    /*if (argc != 11) {
        validInput = 0;
    }*/
    ll SEED = -1;
    ll GIRTH = -1;
    vector<vector<int> > PROTOGRAPH;
    ll CIRCULANT_SIZE = -1;
    ll VARIABLE_NODES;
    ll CHECK_NODES;
    ll DESIRED_NUMBER_OF_MATRICES = -1;
    string INPUT_FILENAME = "";
    for (int i = 1; i + 1 < argc; ++i) {
        if (string(argv[i]) == "-seed") {
            validInput = validInput && toUnsignedInt(argv[i + 1], SEED);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-girth") {
            validInput = validInput && toUnsignedInt(argv[i + 1], GIRTH);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-circulant") {
            validInput = validInput && toUnsignedInt(argv[i + 1], CIRCULANT_SIZE);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-numberOfMatrices") {
            validInput = validInput && toUnsignedInt(argv[i + 1], DESIRED_NUMBER_OF_MATRICES);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-file") {
            INPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-regular") {
            validInput = validInput && toUnsignedInt(argv[i + 1], VARIABLE_NODES) && (i + 2 < argc) && toUnsignedInt(argv[i + 2], CHECK_NODES);
            i += 2;
            regular = 1;
            continue;
        }


    }
    if ((GIRTH < 0) || (SEED < 0) || (CIRCULANT_SIZE < 0) || (DESIRED_NUMBER_OF_MATRICES < 0) || ((INPUT_FILENAME == "") && (!regular)))
        validInput = 0;
    if (!validInput) {
        std::cerr << "Usage: " << argv[0] << " -seed SEED -girth GIRTH -circulant CIRCULANT_SIZE -numberOfMatrices DESIRED_NUMBER_OF_MATRICES -file INPUT_FILENAME" << std::endl;
        return 1;
    }
    if (GIRTH & 1) {
        cerr << "girth must be even\n";
        return 1;
    }
    srand(SEED);
    if (!regular) {
        freopen(INPUT_FILENAME.c_str(), "r", stdin);
        cin >> VARIABLE_NODES >> CHECK_NODES;
        PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES));
        for (int i = 0; i < CHECK_NODES; ++i) {
            for (int j = 0; j < VARIABLE_NODES; ++j) {
                cin >> PROTOGRAPH[i][j];
            }
        }
        fclose(stdin);
    }
    else {
        PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES, 1));
    }

    string folderName = toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH);
    string outputFilename = folderName + "/" + toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH) + "seed" + toStr(SEED);
    if (regular)
        outputFilename += "regular_protograph_matrix"; 
    else
        outputFilename += "protograph_from_" + INPUT_FILENAME + "_matrix";
    if (!isPossible(GIRTH, PROTOGRAPH, CIRCULANT_SIZE))
        return 0;
    system(("mkdir " + folderName).c_str());
    time_t start = time(NULL);
    ll iterationCount = 0;
    ll successCount = 0;
    ll power = 1;
    while (successCount < DESIRED_NUMBER_OF_MATRICES) {
        ++iterationCount;
        vector<vector<vector<int> > > a;
        if (gen(CHECK_NODES, VARIABLE_NODES, PROTOGRAPH, CIRCULANT_SIZE, GIRTH, a)) {
            freopen((outputFilename + toStr(iterationCount) + ".txt").c_str(), "w", stdout);
            ++successCount;
            cout << VARIABLE_NODES << "\t" << CHECK_NODES << "\t" << CIRCULANT_SIZE << endl;
            print(a);
            cout << endl;
            fclose(stdout);
            cerr << "girth = " << getGirth(a, PROTOGRAPH, CIRCULANT_SIZE) << endl;
            eprint(a);
            cerr << iterationCount << " iterations\n";
            cerr << time(NULL) - start << " seconds\n";
            cerr << endl;
        }
        if (iterationCount == power) {
            power *= 2;
            cerr << iterationCount << " iterations\n";
            cerr << time(NULL) - start << " seconds\n";
            cerr << endl;
        }
    }
    return 0;
}