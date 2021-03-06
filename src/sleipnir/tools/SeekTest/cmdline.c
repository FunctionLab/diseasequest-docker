/*
  File autogenerated by gengetopt version 2.22.5
  generated with the following command:
  /usr/bin/gengetopt -iSeekTest.ggo --default-optional -u -N -e 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "Statistical test on a gene set for a given dataset";

const char *gengetopt_args_info_usage = "Usage: SeekTest [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "      --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "\nMode:",
  "  -D, --dab                     DAB mode  (default=off)",
  "  -A, --bin                     PCL Bin mode  (default=off)",
  "  -d, --db                      DB mode  (default=off)",
  "\nDB mode:",
  "  -E, --db_dir=directory        DB directory",
  "  -b, --db_num=INT              Number of files in DB directory  \n                                  (default=`1000')",
  "  -P, --prep=filename           Prep directory (containing .gavg and .gpres \n                                  files)",
  "  -s, --sinfo=filename          Sinfo directory (containing .sinfo files)",
  "  -C, --dataset_list=filename   Dataset list",
  "  -Q, --query=filename          List of genes separated by spaces in one line",
  "  -q, --quant=filename          Quant file",
  "  -c, --count_pair              Count number of z-scores exceeding a threshold  \n                                  (default=off)",
  "  -h, --histogram               Get distribution of z-scores of given genes  \n                                  (default=off)",
  "\nDAB mode:",
  "  -g, --gene_set_list=filename  List of gene-set files",
  "  -x, --input=filename          Gene mapping file",
  "  -B, --dabinput=filename       DAB dataset file",
  "  -a, --gavg_input=filename     Gene average (.gavg) input file",
  "  -p, --gpres_input=filename    Gene presence (.gpres) input file",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->dab_given = 0 ;
  args_info->bin_given = 0 ;
  args_info->db_given = 0 ;
  args_info->db_dir_given = 0 ;
  args_info->db_num_given = 0 ;
  args_info->prep_given = 0 ;
  args_info->sinfo_given = 0 ;
  args_info->dataset_list_given = 0 ;
  args_info->query_given = 0 ;
  args_info->quant_given = 0 ;
  args_info->count_pair_given = 0 ;
  args_info->histogram_given = 0 ;
  args_info->gene_set_list_given = 0 ;
  args_info->input_given = 0 ;
  args_info->dabinput_given = 0 ;
  args_info->gavg_input_given = 0 ;
  args_info->gpres_input_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->dab_flag = 0;
  args_info->bin_flag = 0;
  args_info->db_flag = 0;
  args_info->db_dir_arg = NULL;
  args_info->db_dir_orig = NULL;
  args_info->db_num_arg = 1000;
  args_info->db_num_orig = NULL;
  args_info->prep_arg = NULL;
  args_info->prep_orig = NULL;
  args_info->sinfo_arg = NULL;
  args_info->sinfo_orig = NULL;
  args_info->dataset_list_arg = NULL;
  args_info->dataset_list_orig = NULL;
  args_info->query_arg = NULL;
  args_info->query_orig = NULL;
  args_info->quant_arg = NULL;
  args_info->quant_orig = NULL;
  args_info->count_pair_flag = 0;
  args_info->histogram_flag = 0;
  args_info->gene_set_list_arg = NULL;
  args_info->gene_set_list_orig = NULL;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->dabinput_arg = NULL;
  args_info->dabinput_orig = NULL;
  args_info->gavg_input_arg = NULL;
  args_info->gavg_input_orig = NULL;
  args_info->gpres_input_arg = NULL;
  args_info->gpres_input_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->dab_help = gengetopt_args_info_help[3] ;
  args_info->bin_help = gengetopt_args_info_help[4] ;
  args_info->db_help = gengetopt_args_info_help[5] ;
  args_info->db_dir_help = gengetopt_args_info_help[7] ;
  args_info->db_num_help = gengetopt_args_info_help[8] ;
  args_info->prep_help = gengetopt_args_info_help[9] ;
  args_info->sinfo_help = gengetopt_args_info_help[10] ;
  args_info->dataset_list_help = gengetopt_args_info_help[11] ;
  args_info->query_help = gengetopt_args_info_help[12] ;
  args_info->quant_help = gengetopt_args_info_help[13] ;
  args_info->count_pair_help = gengetopt_args_info_help[14] ;
  args_info->histogram_help = gengetopt_args_info_help[15] ;
  args_info->gene_set_list_help = gengetopt_args_info_help[17] ;
  args_info->input_help = gengetopt_args_info_help[18] ;
  args_info->dabinput_help = gengetopt_args_info_help[19] ;
  args_info->gavg_input_help = gengetopt_args_info_help[20] ;
  args_info->gpres_input_help = gengetopt_args_info_help[21] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->db_dir_arg));
  free_string_field (&(args_info->db_dir_orig));
  free_string_field (&(args_info->db_num_orig));
  free_string_field (&(args_info->prep_arg));
  free_string_field (&(args_info->prep_orig));
  free_string_field (&(args_info->sinfo_arg));
  free_string_field (&(args_info->sinfo_orig));
  free_string_field (&(args_info->dataset_list_arg));
  free_string_field (&(args_info->dataset_list_orig));
  free_string_field (&(args_info->query_arg));
  free_string_field (&(args_info->query_orig));
  free_string_field (&(args_info->quant_arg));
  free_string_field (&(args_info->quant_orig));
  free_string_field (&(args_info->gene_set_list_arg));
  free_string_field (&(args_info->gene_set_list_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->dabinput_arg));
  free_string_field (&(args_info->dabinput_orig));
  free_string_field (&(args_info->gavg_input_arg));
  free_string_field (&(args_info->gavg_input_orig));
  free_string_field (&(args_info->gpres_input_arg));
  free_string_field (&(args_info->gpres_input_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->dab_given)
    write_into_file(outfile, "dab", 0, 0 );
  if (args_info->bin_given)
    write_into_file(outfile, "bin", 0, 0 );
  if (args_info->db_given)
    write_into_file(outfile, "db", 0, 0 );
  if (args_info->db_dir_given)
    write_into_file(outfile, "db_dir", args_info->db_dir_orig, 0);
  if (args_info->db_num_given)
    write_into_file(outfile, "db_num", args_info->db_num_orig, 0);
  if (args_info->prep_given)
    write_into_file(outfile, "prep", args_info->prep_orig, 0);
  if (args_info->sinfo_given)
    write_into_file(outfile, "sinfo", args_info->sinfo_orig, 0);
  if (args_info->dataset_list_given)
    write_into_file(outfile, "dataset_list", args_info->dataset_list_orig, 0);
  if (args_info->query_given)
    write_into_file(outfile, "query", args_info->query_orig, 0);
  if (args_info->quant_given)
    write_into_file(outfile, "quant", args_info->quant_orig, 0);
  if (args_info->count_pair_given)
    write_into_file(outfile, "count_pair", 0, 0 );
  if (args_info->histogram_given)
    write_into_file(outfile, "histogram", 0, 0 );
  if (args_info->gene_set_list_given)
    write_into_file(outfile, "gene_set_list", args_info->gene_set_list_orig, 0);
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->dabinput_given)
    write_into_file(outfile, "dabinput", args_info->dabinput_orig, 0);
  if (args_info->gavg_input_given)
    write_into_file(outfile, "gavg_input", args_info->gavg_input_orig, 0);
  if (args_info->gpres_input_given)
    write_into_file(outfile, "gpres_input", args_info->gpres_input_orig, 0);
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char **argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char **argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  return result;
}

int
cmdline_parser2 (int argc, char **argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (
  int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 0 },
        { "version",	0, NULL, 'V' },
        { "dab",	0, NULL, 'D' },
        { "bin",	0, NULL, 'A' },
        { "db",	0, NULL, 'd' },
        { "db_dir",	1, NULL, 'E' },
        { "db_num",	1, NULL, 'b' },
        { "prep",	1, NULL, 'P' },
        { "sinfo",	1, NULL, 's' },
        { "dataset_list",	1, NULL, 'C' },
        { "query",	1, NULL, 'Q' },
        { "quant",	1, NULL, 'q' },
        { "count_pair",	0, NULL, 'c' },
        { "histogram",	0, NULL, 'h' },
        { "gene_set_list",	1, NULL, 'g' },
        { "input",	1, NULL, 'x' },
        { "dabinput",	1, NULL, 'B' },
        { "gavg_input",	1, NULL, 'a' },
        { "gpres_input",	1, NULL, 'p' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "VDAdE:b:P:s:C:Q:q:chg:x:B:a:p:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'V':	/* Print version and exit.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->version_given),
              &(local_args_info.version_given), optarg, 0, 0, ARG_NO,
              check_ambiguity, override, 0, 0,
              "version", 'V',
              additional_error))
            goto failure;
          cmdline_parser_free (&local_args_info);
          return 0;
        
          break;
        case 'D':	/* DAB mode.  */
        
        
          if (update_arg((void *)&(args_info->dab_flag), 0, &(args_info->dab_given),
              &(local_args_info.dab_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "dab", 'D',
              additional_error))
            goto failure;
        
          break;
        case 'A':	/* PCL Bin mode.  */
        
        
          if (update_arg((void *)&(args_info->bin_flag), 0, &(args_info->bin_given),
              &(local_args_info.bin_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "bin", 'A',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* DB mode.  */
        
        
          if (update_arg((void *)&(args_info->db_flag), 0, &(args_info->db_given),
              &(local_args_info.db_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "db", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'E':	/* DB directory.  */
        
        
          if (update_arg( (void *)&(args_info->db_dir_arg), 
               &(args_info->db_dir_orig), &(args_info->db_dir_given),
              &(local_args_info.db_dir_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "db_dir", 'E',
              additional_error))
            goto failure;
        
          break;
        case 'b':	/* Number of files in DB directory.  */
        
        
          if (update_arg( (void *)&(args_info->db_num_arg), 
               &(args_info->db_num_orig), &(args_info->db_num_given),
              &(local_args_info.db_num_given), optarg, 0, "1000", ARG_INT,
              check_ambiguity, override, 0, 0,
              "db_num", 'b',
              additional_error))
            goto failure;
        
          break;
        case 'P':	/* Prep directory (containing .gavg and .gpres files).  */
        
        
          if (update_arg( (void *)&(args_info->prep_arg), 
               &(args_info->prep_orig), &(args_info->prep_given),
              &(local_args_info.prep_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "prep", 'P',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Sinfo directory (containing .sinfo files).  */
        
        
          if (update_arg( (void *)&(args_info->sinfo_arg), 
               &(args_info->sinfo_orig), &(args_info->sinfo_given),
              &(local_args_info.sinfo_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "sinfo", 's',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Dataset list.  */
        
        
          if (update_arg( (void *)&(args_info->dataset_list_arg), 
               &(args_info->dataset_list_orig), &(args_info->dataset_list_given),
              &(local_args_info.dataset_list_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dataset_list", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'Q':	/* List of genes separated by spaces in one line.  */
        
        
          if (update_arg( (void *)&(args_info->query_arg), 
               &(args_info->query_orig), &(args_info->query_given),
              &(local_args_info.query_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "query", 'Q',
              additional_error))
            goto failure;
        
          break;
        case 'q':	/* Quant file.  */
        
        
          if (update_arg( (void *)&(args_info->quant_arg), 
               &(args_info->quant_orig), &(args_info->quant_given),
              &(local_args_info.quant_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "quant", 'q',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Count number of z-scores exceeding a threshold.  */
        
        
          if (update_arg((void *)&(args_info->count_pair_flag), 0, &(args_info->count_pair_given),
              &(local_args_info.count_pair_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "count_pair", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'h':	/* Get distribution of z-scores of given genes.  */
        
        
          if (update_arg((void *)&(args_info->histogram_flag), 0, &(args_info->histogram_given),
              &(local_args_info.histogram_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "histogram", 'h',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* List of gene-set files.  */
        
        
          if (update_arg( (void *)&(args_info->gene_set_list_arg), 
               &(args_info->gene_set_list_orig), &(args_info->gene_set_list_given),
              &(local_args_info.gene_set_list_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "gene_set_list", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'x':	/* Gene mapping file.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'x',
              additional_error))
            goto failure;
        
          break;
        case 'B':	/* DAB dataset file.  */
        
        
          if (update_arg( (void *)&(args_info->dabinput_arg), 
               &(args_info->dabinput_orig), &(args_info->dabinput_given),
              &(local_args_info.dabinput_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dabinput", 'B',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* Gene average (.gavg) input file.  */
        
        
          if (update_arg( (void *)&(args_info->gavg_input_arg), 
               &(args_info->gavg_input_orig), &(args_info->gavg_input_given),
              &(local_args_info.gavg_input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "gavg_input", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Gene presence (.gpres) input file.  */
        
        
          if (update_arg( (void *)&(args_info->gpres_input_arg), 
               &(args_info->gpres_input_orig), &(args_info->gpres_input_given),
              &(local_args_info.gpres_input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "gpres_input", 'p',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          if (strcmp (long_options[option_index].name, "help") == 0) {
            cmdline_parser_print_help ();
            cmdline_parser_free (&local_args_info);
            exit (EXIT_SUCCESS);
          }

        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */




  cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
