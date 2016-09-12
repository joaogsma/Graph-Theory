1. O arquivo main.cpp usa funcionalidades de C++11, então para compilá-lo, talvez seja necessário habilitar alguma flag no compilador, já que a maioria dos compiladores requer isso. Na minha máquina (Windows), eu utilizei o compilador g++ e foi necessário habilitar a flag -std=c++0x. O comando de compilação inteiro era "g++ ProjetoMF_jgsma.cpp -std=c++0x".

2. Não foi especificado como os números nos dados de entrada deveriam ser separados, então eu assumi que são separados por espaços. Um exemplo de conteúdo de um arquivo de entrada:
1
10 7
1 2 3
1 4 3
2 3 4
3 4 1
3 5 2
5 2 1
4 5 2
4 6 6
5 7 1
6 7 9