
    #include <iostream>
    #include <fstream>
    #include <vector>
    #include <numeric>
    #include <queue>
    #include <algorithm>
    #include <functional>
    #include <iomanip>
    #include <string>
    #include <stack>
    #include <random>
    #include <chrono>
    using namespace std;

    enum Representacao { LISTA, MATRIZ };

    struct Node {
        int v;
        Node* next;
        Node(int vert, Node* nxt = nullptr) : v(vert), next(nxt) {}
    };

    class Grafo {
    private:
        int n; // número de vértices
        int m; // número de arestas
        Representacao tipo;
        vector<Node*> adjList;        // lista encadeada
        vector<vector<int>> adjMatrix; // matriz de adjacência

    public:
        Grafo(int x, Representacao r) : n(x), m(0), tipo(r) {
            if (tipo == LISTA) {
                adjList.assign(n + 1, nullptr);
            } else {
                adjMatrix.assign(n + 1, vector<int>(n + 1, 0));
            }
        }

        void adicionarAresta(int u, int v) {
            if (tipo == LISTA) {
                adjList[u] = new Node(v, adjList[u]);
                adjList[v] = new Node(u, adjList[v]);
            } else {
                adjMatrix[u][v] = adjMatrix[v][u] = 1;
            }
            m++;
        }

        static Grafo lerDeArquivo(const string& nomeArquivo, Representacao r) {
            ifstream entrada(nomeArquivo);
            if (!entrada.is_open()) {
                cerr << "Erro ao abrir o arquivo." << endl;
                exit(1);
            }
            int vertices, u, v;
            entrada >> vertices;
            Grafo g(vertices, r);
            while (entrada >> u >> v) {
                g.adicionarAresta(u, v);
            }
            entrada.close();
            return g;
        }

        int contarArestas() const {
            if (tipo == LISTA) {
                int total = 0;
                for (int i = 1; i <= n; i++) {
                    Node* curr = adjList[i];
                    while (curr) {
                        total++;
                        curr = curr->next;
                    }
                }
                return total / 2; // cada aresta foi contada 2 vezes
            } else {
                int total = 0;
                for (int i = 1; i <= n; i++)
                    for (int j = i + 1; j <= n; j++)
                        if (adjMatrix[i][j]) total++;
                return total;
            }
        }

        vector<int> graus() const {
            vector<int> g(n + 1, 0);
            if (tipo == LISTA) {
                for (int i = 1; i <= n; i++) {
                    int grau = 0;
                    Node* curr = adjList[i];
                    while (curr) {
                        grau++;
                        curr = curr->next;
                    }
                    g[i] = grau;
                }
            } else {
                for (int i = 1; i <= n; i++) {
                    int grau = 0;
                    for (int j = 1; j <= n; j++) if (adjMatrix[i][j]) grau++;
                    g[i] = grau;
                }
            }
            return g;
        }

        int distancia(int origem, int destino) const {
            if (origem < 1 || origem > n || destino < 1 || destino > n) return -1;
            vector<int> dist(n + 1, -1);
            queue<int> q;
            dist[origem] = 0;
            q.push(origem);
            while (!q.empty()) {
                int u = q.front(); q.pop();
                if (tipo == LISTA) {
                    for (Node* curr = adjList[u]; curr; curr = curr->next) {
                        int v = curr->v;
                        if (dist[v] == -1) {
                            dist[v] = dist[u] + 1;
                            if (v == destino) return dist[v];
                            q.push(v);
                        }
                    }
                } else {
                    for (int v = 1; v <= n; v++) {
                        if (adjMatrix[u][v] && dist[v] == -1) {
                            dist[v] = dist[u] + 1;
                            if (v == destino) return dist[v];
                            q.push(v);
                        }
                    }
                }
            }
            return -1;
        }

        int diametro() const {
            int diam = 0;
            for (int i = 1; i <= n; i++) {
                vector<int> dist(n + 1, -1);
                queue<int> q;
                dist[i] = 0;
                q.push(i);
                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    if (tipo == LISTA) {
                        for (Node* curr = adjList[u]; curr; curr = curr->next) {
                            int v = curr->v;
                            if (dist[v] == -1) {
                                dist[v] = dist[u] + 1;
                                diam = max(diam, dist[v]);
                                q.push(v);
                            }
                        }
                    } else {
                        for (int v = 1; v <= n; v++) {
                            if (adjMatrix[u][v] && dist[v] == -1) {
                                dist[v] = dist[u] + 1;
                                diam = max(diam, dist[v]);
                                q.push(v);
                            }
                        }
                    }
                }
            }
            return diam;
        }

        vector<vector<int>> componentesConexas() const {
            vector<bool> visitado(n + 1, false);
            vector<vector<int>> componentes;
            for (int i = 1; i <= n; i++) {
                if (!visitado[i]) {
                    vector<int> componente;
                    queue<int> q;
                    q.push(i);
                    visitado[i] = true;
                    while (!q.empty()) {
                        int u = q.front(); q.pop();
                        componente.push_back(u);
                        if (tipo == LISTA) {
                            for (Node* curr = adjList[u]; curr; curr = curr->next) {
                                int v = curr->v;
                                if (!visitado[v]) { visitado[v] = true; q.push(v); }
                            }
                        } else {
                            for (int v = 1; v <= n; v++) {
                                if (adjMatrix[u][v] && !visitado[v]) { visitado[v] = true; q.push(v); }
                            }
                        }
                    }
                    componentes.push_back(componente);
                }
            }
            sort(componentes.begin(), componentes.end(),
                [](const vector<int>& a, const vector<int>& b){ return a.size() > b.size(); });
            return componentes;
        }

        void salvarEstatisticas(const string& nomeArquivo) const {
            ofstream out(nomeArquivo);
            if (!out.is_open()) {
                throw runtime_error("Erro ao criar arquivo de saída!");
            }

            int nArestas = contarArestas();
            vector<int> gAll = graus();
            vector<int> g = vector<int>(gAll.begin() + 1, gAll.end());

            int grauMin = *min_element(g.begin(), g.end());
            int grauMax = *max_element(g.begin(), g.end());
            double grauMedio = accumulate(g.begin(), g.end(), 0.0) / (double)n;

            vector<int> ordenado = g;
            sort(ordenado.begin(), ordenado.end());
            double mediana;
            if (n % 2 == 0) mediana = (ordenado[n/2 - 1] + ordenado[n/2]) / 2.0;
            else mediana = ordenado[n/2];

            int diam = diametro();
            auto comps = componentesConexas();

            out << fixed << setprecision(2);
            out << "Numero de vertices: " << n << "\n";
            out << "Numero de arestas: " << nArestas << "\n";
            out << "Grau minimo: " << grauMin << "\n";
            out << "Grau maximo: " << grauMax << "\n";
            out << "Grau medio: " << grauMedio << "\n";
            out << "Mediana do grau: " << mediana << "\n";
            out << "Diametro do grafo: " << diam << "\n\n";

            out << "Numero de componentes conexas: " << comps.size() << "\n\n";
            for (size_t i = 0; i < comps.size(); ++i) {
                out << "Componente " << i+1 << " (tamanho " << comps[i].size() << "): ";
                for (int v : comps[i]) out << v << " ";
                out << "\n";
            }

            out.close();
        }

        int diametroAproximado() const {
        auto bfsMaisDistante = [&](int origem) {
            vector<int> dist(n + 1, -1);
            queue<int> q;
            dist[origem] = 0;
            q.push(origem);
            int ultimo = origem;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                ultimo = u;
                if (tipo == LISTA) {
                    for (Node* curr = adjList[u]; curr; curr = curr->next) {
                        int v = curr->v;
                        if (dist[v] == -1) {
                            dist[v] = dist[u] + 1;
                            q.push(v);
                        }
                    }
                } else {
                    for (int v = 1; v <= n; v++) {
                        if (adjMatrix[u][v] && dist[v] == -1) {
                            dist[v] = dist[u] + 1;
                            q.push(v);
                        }
                    }
                }
            }
            return pair<int,int>(ultimo, dist[ultimo]);
        };

        int inicio = 1;
        pair<int, int> res1 = bfsMaisDistante(inicio);
        pair<int, int> res2 = bfsMaisDistante(res1.first);
        return res2.second;
    }

        void BFS(int origem, const string& arquivo) const {
            vector<int> nivel(n + 1, -1);
            vector<int> pai(n + 1, 0);
            queue<int> q;

            nivel[origem] = 0;
            pai[origem] = 0;
            q.push(origem);

            while (!q.empty()) {
                int u = q.front(); q.pop();
                if (tipo == LISTA) {
                    for (Node* curr = adjList[u]; curr; curr = curr->next) {
                        int v = curr->v;
                        if (nivel[v] == -1) {
                            pai[v] = u;
                            nivel[v] = nivel[u] + 1;
                            q.push(v);
                        }
                    }
                } else {
                    for (int v = 1; v <= n; v++) {
                        if (adjMatrix[u][v] && nivel[v] == -1) {
                            pai[v] = u;
                            nivel[v] = nivel[u] + 1;
                            q.push(v);
                        }
                    }
                }
            }

            ofstream out(arquivo);
            out << "BFS - Arvore de busca a partir do vertice " << origem << "\n";
            out << "Vertice\tPai\tNivel\n";
            for (int i = 1; i <= n; i++) {
                out << i << "\t\t" << pai[i] << "\t" << nivel[i] << "\n";
            }
            out.close();
        }

        void DFS(int origem, const string& arquivo) const {
            vector<int> nivel(n + 1, -1);
            vector<int> pai(n + 1, 0);
            stack<pair<int,int>> st;

            st.push({origem, 0});
            pai[origem] = 0;

            while (!st.empty()) {
                auto [u, d] = st.top();
                st.pop();

                if (nivel[u] != -1) continue; // já visitado
                nivel[u] = d;

                if (tipo == LISTA) {
                    vector<int> vizinhos;
                    for (Node* curr = adjList[u]; curr; curr = curr->next) {
                        vizinhos.push_back(curr->v);
                    }
                    sort(vizinhos.rbegin(), vizinhos.rend()); 
                    for (int v : vizinhos) {
                        if (nivel[v] == -1) {
                            pai[v] = u;
                            st.push({v, d + 1});
                        }
                    }
                } else { // MATRIZ
                    for (int v = n; v >= 1; v--) {
                        if (adjMatrix[u][v] && nivel[v] == -1) {
                            pai[v] = u;
                            st.push({v, d + 1});
                        }
                    }
                }
            }

            ofstream out(arquivo);
            out << "DFS - Arvore de busca a partir do vertice " << origem << "\n";
            out << "Vertice\tPai\tNivel\n";
            for (int i = 1; i <= n; i++) {
                out << i << "\t\t" << pai[i] << "\t" << nivel[i] << "\n";
            }
            out.close();
        }

        ~Grafo() {
            if (tipo == LISTA) {
                for (int i = 1; i <= n; i++) {
                    Node* curr = adjList[i];
                    while (curr) {
                        Node* tmp = curr;
                        curr = curr->next;
                        delete tmp;
                    }
                }
            }
        }
    };

int main() {
    try {
        Grafo gLista = Grafo::lerDeArquivo("teste.txt", LISTA);

        // === Caso 6 ===
        auto componentes = gLista.componentesConexas();
        cout << "Numero de componentes conexas: " << componentes.size() << "\n";
        if (!componentes.empty()) {
            cout << "Tamanho da maior componente: " << componentes.front().size() << "\n";
            cout << "Tamanho da menor componente: " << componentes.back().size() << "\n";
        }

        // === Caso 7 ===
        int diam = gLista.diametro();
        int diamApprox = gLista.diametroAproximado();
        cout << "Diametro exato: " << diam << "\n";
        cout << "Diametro aproximado: " << diamApprox << "\n";

    } catch (const exception& e) {
        cerr << "Erro: " << e.what() << "\n";
    }
    return 0;
}