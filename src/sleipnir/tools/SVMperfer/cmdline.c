/*
  File autogenerated by gengetopt version 2.22.5
  generated with the following command:
  /usr/bin/gengetopt -iSVMperfer.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Wrapper for SVM perf";

const char *gengetopt_args_info_usage = "Usage: SVMperfer [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                  Print help and exit",
  "  -V, --version               Print version and exit",
  "\nMain:",
  "  -l, --labels=filename       Labels file",
  "  -o, --output=filename       Output file ",
  "  -i, --input=filename        Input PCL file ",
  "  -m, --model=filename        Model file",
  "  -T, --test_labels=filename  Test Labels file",
  "  -a, --all                   Always classify all genes in PCLs  (default=off)",
  "  -S, --slack                 Use slack rescaling (not implemented for ROC \n                                loss)  (default=off)",
  "\nOptions:",
  "  -v, --verbosity=INT         Sets the svm_struct verbosity  (default=`0')",
  "  -s, --skip=INT              Number of columns to skip in input pcls  \n                                (default=`2')",
  "  -n, --normalize             Normalize PCLS to 0 mean 1 variance  \n                                (default=off)",
  "  -c, --cross_validation=INT  Number of cross-validation sets ( arg of 1 will \n                                turn off cross-validation )  (default=`5')",
  "  -e, --error_function=INT    Sets the loss function for SVM learning: Choice \n                                of:\n\n                                0\tZero/one loss: 1 if vector of predictions \n                                contains error, 0 otherwise.\n\n                                1\tF1: 100 minus the F1-score in percent.\n\n                                2\tErrorrate: Percentage of errors in \n                                prediction vector.\n\n                                3\tPrec/Rec Breakeven: 100 minus PRBEP in \n                                percent.\n\n                                4\tPrec@k: 100 minus precision at k in percent.\n\n                                5\tRec@k: 100 minus recall at k in percent.\n\n                                10\tROCArea: Percentage of swapped pos/neg \n                                pairs (i.e. 100 - ROCArea).\n                                  (default=`10')",
  "  -k, --k_value=FLOAT         Value of k parameter used for Prec@k and Rec@k in \n                                (0,1)  (default=`0.5')",
  "  -t, --tradeoff=FLOAT        SVM tradeoff constant C  (default=`1')",
  "  -A, --simple_model          Write model files with only linear weights  \n                                (default=on)",
  "  -p, --params=filename       Parameter file",
  "  -M, --mmap                  Memory map binary input  (default=off)",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_FLOAT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);

static int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->labels_given = 0 ;
  args_info->output_given = 0 ;
  args_info->input_given = 0 ;
  args_info->model_given = 0 ;
  args_info->test_labels_given = 0 ;
  args_info->all_given = 0 ;
  args_info->slack_given = 0 ;
  args_info->verbosity_given = 0 ;
  args_info->skip_given = 0 ;
  args_info->normalize_given = 0 ;
  args_info->cross_validation_given = 0 ;
  args_info->error_function_given = 0 ;
  args_info->k_value_given = 0 ;
  args_info->tradeoff_given = 0 ;
  args_info->simple_model_given = 0 ;
  args_info->params_given = 0 ;
  args_info->mmap_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->labels_arg = NULL;
  args_info->labels_orig = NULL;
  args_info->output_arg = NULL;
  args_info->output_orig = NULL;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->model_arg = NULL;
  args_info->model_orig = NULL;
  args_info->test_labels_arg = NULL;
  args_info->test_labels_orig = NULL;
  args_info->all_flag = 0;
  args_info->slack_flag = 0;
  args_info->verbosity_arg = 0;
  args_info->verbosity_orig = NULL;
  args_info->skip_arg = 2;
  args_info->skip_orig = NULL;
  args_info->normalize_flag = 0;
  args_info->cross_validation_arg = 5;
  args_info->cross_validation_orig = NULL;
  args_info->error_function_arg = 10;
  args_info->error_function_orig = NULL;
  args_info->k_value_arg = 0.5;
  args_info->k_value_orig = NULL;
  args_info->tradeoff_arg = 1;
  args_info->tradeoff_orig = NULL;
  args_info->simple_model_flag = 1;
  args_info->params_arg = NULL;
  args_info->params_orig = NULL;
  args_info->mmap_flag = 0;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->labels_help = gengetopt_args_info_help[3] ;
  args_info->output_help = gengetopt_args_info_help[4] ;
  args_info->input_help = gengetopt_args_info_help[5] ;
  args_info->model_help = gengetopt_args_info_help[6] ;
  args_info->test_labels_help = gengetopt_args_info_help[7] ;
  args_info->all_help = gengetopt_args_info_help[8] ;
  args_info->slack_help = gengetopt_args_info_help[9] ;
  args_info->verbosity_help = gengetopt_args_info_help[11] ;
  args_info->skip_help = gengetopt_args_info_help[12] ;
  args_info->normalize_help = gengetopt_args_info_help[13] ;
  args_info->cross_validation_help = gengetopt_args_info_help[14] ;
  args_info->error_function_help = gengetopt_args_info_help[15] ;
  args_info->k_value_help = gengetopt_args_info_help[16] ;
  args_info->tradeoff_help = gengetopt_args_info_help[17] ;
  args_info->simple_model_help = gengetopt_args_info_help[18] ;
  args_info->params_help = gengetopt_args_info_help[19] ;
  args_info->mmap_help = gengetopt_args_info_help[20] ;
  
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
  free_string_field (&(args_info->labels_arg));
  free_string_field (&(args_info->labels_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->model_arg));
  free_string_field (&(args_info->model_orig));
  free_string_field (&(args_info->test_labels_arg));
  free_string_field (&(args_info->test_labels_orig));
  free_string_field (&(args_info->verbosity_orig));
  free_string_field (&(args_info->skip_orig));
  free_string_field (&(args_info->cross_validation_orig));
  free_string_field (&(args_info->error_function_orig));
  free_string_field (&(args_info->k_value_orig));
  free_string_field (&(args_info->tradeoff_orig));
  free_string_field (&(args_info->params_arg));
  free_string_field (&(args_info->params_orig));
  
  
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
  if (args_info->labels_given)
    write_into_file(outfile, "labels", args_info->labels_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->model_given)
    write_into_file(outfile, "model", args_info->model_orig, 0);
  if (args_info->test_labels_given)
    write_into_file(outfile, "test_labels", args_info->test_labels_orig, 0);
  if (args_info->all_given)
    write_into_file(outfile, "all", 0, 0 );
  if (args_info->slack_given)
    write_into_file(outfile, "slack", 0, 0 );
  if (args_info->verbosity_given)
    write_into_file(outfile, "verbosity", args_info->verbosity_orig, 0);
  if (args_info->skip_given)
    write_into_file(outfile, "skip", args_info->skip_orig, 0);
  if (args_info->normalize_given)
    write_into_file(outfile, "normalize", 0, 0 );
  if (args_info->cross_validation_given)
    write_into_file(outfile, "cross_validation", args_info->cross_validation_orig, 0);
  if (args_info->error_function_given)
    write_into_file(outfile, "error_function", args_info->error_function_orig, 0);
  if (args_info->k_value_given)
    write_into_file(outfile, "k_value", args_info->k_value_orig, 0);
  if (args_info->tradeoff_given)
    write_into_file(outfile, "tradeoff", args_info->tradeoff_orig, 0);
  if (args_info->simple_model_given)
    write_into_file(outfile, "simple_model", 0, 0 );
  if (args_info->params_given)
    write_into_file(outfile, "params", args_info->params_orig, 0);
  if (args_info->mmap_given)
    write_into_file(outfile, "mmap", 0, 0 );
  

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
  int result = EXIT_SUCCESS;

  if (cmdline_parser_required2(args_info, prog_name, 0) > 0)
    result = EXIT_FAILURE;

  return result;
}

int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;
  FIX_UNUSED (additional_error);

  /* checks for required options */
  if (! args_info->input_given)
    {
      fprintf (stderr, "%s: '--input' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  
  /* checks for dependences among options */

  return error;
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
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
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
  case ARG_FLOAT:
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
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "labels",	1, NULL, 'l' },
        { "output",	1, NULL, 'o' },
        { "input",	1, NULL, 'i' },
        { "model",	1, NULL, 'm' },
        { "test_labels",	1, NULL, 'T' },
        { "all",	0, NULL, 'a' },
        { "slack",	0, NULL, 'S' },
        { "verbosity",	1, NULL, 'v' },
        { "skip",	1, NULL, 's' },
        { "normalize",	0, NULL, 'n' },
        { "cross_validation",	1, NULL, 'c' },
        { "error_function",	1, NULL, 'e' },
        { "k_value",	1, NULL, 'k' },
        { "tradeoff",	1, NULL, 't' },
        { "simple_model",	0, NULL, 'A' },
        { "params",	1, NULL, 'p' },
        { "mmap",	0, NULL, 'M' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVl:o:i:m:T:aSv:s:nc:e:k:t:Ap:M", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

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
        case 'l':	/* Labels file.  */
        
        
          if (update_arg( (void *)&(args_info->labels_arg), 
               &(args_info->labels_orig), &(args_info->labels_given),
              &(local_args_info.labels_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "labels", 'l',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output file .  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'i':	/* Input PCL file .  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* Model file.  */
        
        
          if (update_arg( (void *)&(args_info->model_arg), 
               &(args_info->model_orig), &(args_info->model_given),
              &(local_args_info.model_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "model", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'T':	/* Test Labels file.  */
        
        
          if (update_arg( (void *)&(args_info->test_labels_arg), 
               &(args_info->test_labels_orig), &(args_info->test_labels_given),
              &(local_args_info.test_labels_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "test_labels", 'T',
              additional_error))
            goto failure;
        
          break;
        case 'a':	/* Always classify all genes in PCLs.  */
        
        
          if (update_arg((void *)&(args_info->all_flag), 0, &(args_info->all_given),
              &(local_args_info.all_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "all", 'a',
              additional_error))
            goto failure;
        
          break;
        case 'S':	/* Use slack rescaling (not implemented for ROC loss).  */
        
        
          if (update_arg((void *)&(args_info->slack_flag), 0, &(args_info->slack_given),
              &(local_args_info.slack_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "slack", 'S',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Sets the svm_struct verbosity.  */
        
        
          if (update_arg( (void *)&(args_info->verbosity_arg), 
               &(args_info->verbosity_orig), &(args_info->verbosity_given),
              &(local_args_info.verbosity_given), optarg, 0, "0", ARG_INT,
              check_ambiguity, override, 0, 0,
              "verbosity", 'v',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Number of columns to skip in input pcls.  */
        
        
          if (update_arg( (void *)&(args_info->skip_arg), 
               &(args_info->skip_orig), &(args_info->skip_given),
              &(local_args_info.skip_given), optarg, 0, "2", ARG_INT,
              check_ambiguity, override, 0, 0,
              "skip", 's',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Normalize PCLS to 0 mean 1 variance.  */
        
        
          if (update_arg((void *)&(args_info->normalize_flag), 0, &(args_info->normalize_given),
              &(local_args_info.normalize_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "normalize", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Number of cross-validation sets ( arg of 1 will turn off cross-validation ).  */
        
        
          if (update_arg( (void *)&(args_info->cross_validation_arg), 
               &(args_info->cross_validation_orig), &(args_info->cross_validation_given),
              &(local_args_info.cross_validation_given), optarg, 0, "5", ARG_INT,
              check_ambiguity, override, 0, 0,
              "cross_validation", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Sets the loss function for SVM learning: Choice of:
        0\tZero/one loss: 1 if vector of predictions contains error, 0 otherwise.
        1\tF1: 100 minus the F1-score in percent.
        2\tErrorrate: Percentage of errors in prediction vector.
        3\tPrec/Rec Breakeven: 100 minus PRBEP in percent.
        4\tPrec@k: 100 minus precision at k in percent.
        5\tRec@k: 100 minus recall at k in percent.
        10\tROCArea: Percentage of swapped pos/neg pairs (i.e. 100 - ROCArea).
.  */
        
        
          if (update_arg( (void *)&(args_info->error_function_arg), 
               &(args_info->error_function_orig), &(args_info->error_function_given),
              &(local_args_info.error_function_given), optarg, 0, "10", ARG_INT,
              check_ambiguity, override, 0, 0,
              "error_function", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'k':	/* Value of k parameter used for Prec@k and Rec@k in (0,1).  */
        
        
          if (update_arg( (void *)&(args_info->k_value_arg), 
               &(args_info->k_value_orig), &(args_info->k_value_given),
              &(local_args_info.k_value_given), optarg, 0, "0.5", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "k_value", 'k',
              additional_error))
            goto failure;
        
          break;
        case 't':	/* SVM tradeoff constant C.  */
        
        
          if (update_arg( (void *)&(args_info->tradeoff_arg), 
               &(args_info->tradeoff_orig), &(args_info->tradeoff_given),
              &(local_args_info.tradeoff_given), optarg, 0, "1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "tradeoff", 't',
              additional_error))
            goto failure;
        
          break;
        case 'A':	/* Write model files with only linear weights.  */
        
        
          if (update_arg((void *)&(args_info->simple_model_flag), 0, &(args_info->simple_model_given),
              &(local_args_info.simple_model_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "simple_model", 'A',
              additional_error))
            goto failure;
        
          break;
        case 'p':	/* Parameter file.  */
        
        
          if (update_arg( (void *)&(args_info->params_arg), 
               &(args_info->params_orig), &(args_info->params_given),
              &(local_args_info.params_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "params", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'M':	/* Memory map binary input.  */
        
        
          if (update_arg((void *)&(args_info->mmap_flag), 0, &(args_info->mmap_given),
              &(local_args_info.mmap_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "mmap", 'M',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */



  if (check_required)
    {
      error += cmdline_parser_required2 (args_info, argv[0], additional_error);
    }

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
