#include <stdio.h>
#include <stdlib.h>

/* 声明 GenoemLib.c 中的函数 */
int RefGene_GetMatchedControl_Main(char strDatabasePath[], int nDatabaseType,
        char strSpecies[], char strChrLenPath[],
        char strInputPath[], char strOutputPath[],
        int nRepNum, int nRegionLen, int nRemoveRedundancy);

int main(int argc, char *argv[])
{
    if (argc < 7) {
        printf("Usage: ./getmatchedcontrol <refgene_db> <species> <chrlen> <input.cod> <output.cod> <regionLen>\n");
        return 1;
    }

    char *refgene_db = argv[1];
    char *species    = argv[2];
    char *chrlen     = argv[3];
    char *input_cod  = argv[4];
    char *output_cod = argv[5];
    int regionLen    = atoi(argv[6]);

    int ret = RefGene_GetMatchedControl_Main(
        refgene_db,      /* refgene database */
        2,               /* nDatabaseType */
        species,         /* hg19 / hg38 / mm10 */
        chrlen,          /* chromosome length file */
        input_cod,       /* peaks.cod */
        output_cod,      /* output cod */
        10,              /* nRepNum */
        regionLen,       /* nRegionLen */
        1                /* remove redundancy */
    );

    printf("Done. Return code = %d\n", ret);
    return 0;
}