tenho que ter as seguintes respostas:
- porque usa intptr_t?
- porque populate virou uma função quer retorna um ponteiro?

algumas mudanças que inicialmente não parecem ser tão importantes:
- var global com o número de threads
- número de threads é inserida pela linha de comando

alguns comentários:
- só paralelizamos o que consideramos principal, deixamos de paralelizar:
    - for que gera o data no main
    - for do kmeans que deixa todos pontos não mapeados
    - for do kmeans que deixa partições sujas e define os seus centros
    - for do kmeans que mapeia pontos não mapeados
- no while do kmeans, cada chamada de *Paralelizado cria threads, da join e libera. Fizemos dessa forma por ser mais fácil e
