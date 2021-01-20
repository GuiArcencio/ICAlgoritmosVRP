# ICAlgoritmosVRP

Repositório de códigos utilizados no projeto de Iniciação Científica que busca atacar o Problema do Roteamento de Veículos com métodos exatos, mais especificamente, Programação Linear Inteira.

## Compilação

Necessário ter o Gurobi instalado. Basta executar:

`make exec`

## Uso

`CVRPSolver [opções...]`

- **-f, --file** \<nome-do-arquivo> | caminho para o arquivo que descreve a instância
- **-H, --use-heuristic** | usa uma solução heurística de \<nome-do-arquivo>.heu
- **C, --karger-coefficient** \<coeficiente> | constante que multiplica o número de vezes que o Karger será executado, 10.0 por padrão
- **-l, --use-log-n** | executa o Karger O(log n) vezes em vez de O(n)
- **-t, --time-limit** \<tempo> | tempo limite para o solver em segundos, 3600.0 por padrão
- **-h, --help** | mostra a página de ajuda

## Plotter

Gera imagens da configuração dos clientes da instância e das rotas geradas. Utilização:

`plot.py <nome-do-arquivo>`

O script utilizará \<nome-do-arquivo> e \<nome-do-arquivo>.sol para a criação das imagens.