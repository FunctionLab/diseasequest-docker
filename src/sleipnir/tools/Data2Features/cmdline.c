/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /home/chuttenh/hg/sleipnir/trunk/../extlib/gengetopt-2.22/bin/gengetopt -iData2Features.ggo --default-optional -u -N -e 

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

#include "getopt.h"

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "Data transformation to feature sets for machine learning";

const char *gengetopt_args_info_usage = "Usage: Data2Features [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                  Print help and exit",
  "  -V, --version               Print version and exit",
  "\nMain:",
  "  -p, --positives=filename    Positive gene list",
  "  -e, --environment=filename  List of environment features and default values",
  "  -d, --data=filename         Feature values for each data set",
  "\nMiscellaneous:",
  "  -g, --genome=filename       SGD features file",
  "\nPCL Processing:",
  "  -D, --distance=STRING       Similarity measure  (possible values=\"pearson\", \n                                \"euclidean\", \"kendalls\", \"kolm-smir\", \n                                \"spearman\", \"pearnorm\", \"hypergeom\", \n                                \"innerprod\", \"bininnerprod\", \"quickpear\", \n                                \"mi\" default=`pearnorm')",
  "  -N, --normalize             Normalize distances  (default=off)",
  "  -Z, --zscore                Convert correlations to z-scores  (default=on)",
  "  -S, --skip=INT              PCL columns to skip after ID  (default=`2')",
  "\nOptional:",
  "  -m, --memmap                Memory map input DABs  (default=off)",
  "  -v, --verbosity=INT         Message verbosity  (default=`5')",
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
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);

static int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);

char *cmdline_parser_distance_values[] = {"pearson", "euclidean", "kendalls", "kolm-smir", "spearman", "pearnorm", "hypergeom", "innerprod", "bininnerprod", "quickpear", "mi", 0} ;	/* Possible values for distance.  */

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->positives_given = 0 ;
  args_info->environment_given = 0 ;
  args_info->data_given = 0 ;
  args_info->genome_given = 0 ;
  args_info->distance_given = 0 ;
  args_info->normalize_given = 0 ;
  args_info->zscore_given = 0 ;
  args_info->skip_given = 0 ;
  args_info->memmap_given = 0 ;
  args_info->verbosity_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->positives_arg = NULL;
  args_info->positives_orig = NULL;
  args_info->environment_arg = NULL;
  args_info->environment_orig = NULL;
  args_info->data_arg = NULL;
  args_info->data_orig = NULL;
  args_info->genome_arg = NULL;
  args_info->genome_orig = NULL;
  args_info->distance_arg = gengetopt_strdup ("pearnorm");
  args_info->distance_orig = NULL;
  args_info->normalize_flag = 0;
  args_info->zscore_flag = 1;
  args_info->skip_arg = 2;
  args_info->skip_orig = NULL;
  args_info->memmap_flag = 0;
  args_info->verbosity_arg = 5;
  args_info->verbosity_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->positives_help = gengetopt_args_info_help[3] ;
  args_info->environment_help = gengetopt_args_info_help[4] ;
  args_info->data_help = gengetopt_args_info_help[5] ;
  args_info->genome_help = gengetopt_args_info_help[7] ;
  args_info->distance_help = gengetopt_args_info_help[9] ;
  args_info->normalize_help = gengetopt_args_info_help[10] ;
  args_info->zscore_help = gengetopt_args_info_help[11] ;
  args_info->skip_help = gengetopt_args_info_help[12] ;
  args_info->memmap_help = gengetopt_args_info_help[14] ;
  args_info->verbosity_help = gengetopt_args_info_help[15] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n", CMDLINE_PARSER_PACKAGE, CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n", gengetopt_args_info_description);
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

  args_info->inputs = NULL;
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
  free_string_field (&(args_info->positives_arg));
  free_string_field (&(args_info->positives_orig));
  free_string_field (&(args_info->environment_arg));
  free_string_field (&(args_info->environment_orig));
  free_string_field (&(args_info->data_arg));
  free_string_field (&(args_info->data_orig));
  free_string_field (&(args_info->genome_arg));
  free_string_field (&(args_info->genome_orig));
  free_string_field (&(args_info->distance_arg));
  free_string_field (&(args_info->distance_orig));
  free_string_field (&(args_info->skip_orig));
  free_string_field (&(args_info->verbosity_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}

/**
 * @param val the value to check
 * @param values the possible values
 * @return the index of the matched value:
 * -1 if no value matched,
 * -2 if more than one value has matched
 */
static int
check_possible_values(const char *val, char *values[])
{
  int i, found, last;
  size_t len;

  if (!val)   /* otherwise strlen() crashes below */
    return -1; /* -1 means no argument for the option */

  found = last = 0;

  for (i = 0, len = strlen(val); values[i]; ++i)
    {
      if (strncmp(val, values[i], len) == 0)
        {
          ++found;
          last = i;
          if (strlen(values[i]) == len)
            return i; /* exact macth no need to check more */
        }
    }

  if (found == 1) /* one match: OK */
    return last;

  return (found ? -2 : -1); /* return many values or none matched */
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, char *values[])
{
  int found = -1;
  if (arg) {
    if (values) {
      found = check_possible_values(arg, values);      
    }
    if (found >= 0)
      fprintf(outfile, "%s=\"%s\" # %s\n", opt, arg, values[found]);
    else
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
  if (args_info->positives_given)
    write_into_file(outfile, "positives", args_info->positives_orig, 0);
  if (args_info->environment_given)
    write_into_file(outfile, "environment", args_info->environment_orig, 0);
  if (args_info->data_given)
    write_into_file(outfile, "data", args_info->data_orig, 0);
  if (args_info->genome_given)
    write_into_file(outfile, "genome", args_info->genome_orig, 0);
  if (args_info->distance_given)
    write_into_file(outfile, "distance", args_info->distance_orig, cmdline_parser_distance_values);
  if (args_info->normalize_given)
    write_into_file(outfile, "normalize", 0, 0 );
  if (args_info->zscore_given)
    write_into_file(outfile, "zscore", 0, 0 );
  if (args_info->skip_given)
    write_into_file(outfile, "skip", args_info->skip_orig, 0);
  if (args_info->memmap_given)
    write_into_file(outfile, "memmap", 0, 0 );
  if (args_info->verbosity_given)
    write_into_file(outfile, "verbosity", args_info->verbosity_orig, 0);
  

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
  char *result = NULL;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char * const *argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char * const *argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, NULL);

  return result;
}

int
cmdline_parser2 (int argc, char * const *argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, NULL);

  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  int result = EXIT_SUCCESS;

  if (cmdline_parser_required2(args_info, prog_name, NULL) > 0)
    result = EXIT_FAILURE;

  return result;
}

int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;

  /* checks for required options */
  if (! args_info->environment_given)
    {
      fprintf (stderr, "%s: '--environment' ('-e') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->data_given)
    {
      fprintf (stderr, "%s: '--data' ('-d') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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
               char *value, char *possible_values[], const char *default_value,
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

  if (possible_values && (found = check_possible_values((value ? value : default_value), possible_values)) < 0)
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: %s argument, \"%s\", for option `--%s' (`-%c')%s\n", 
          package_name, (found == -2) ? "ambiguous" : "invalid", value, long_opt, short_opt,
          (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: %s argument, \"%s\", for option `--%s'%s\n", 
          package_name, (found == -2) ? "ambiguous" : "invalid", value, long_opt,
          (additional_error ? additional_error : ""));
      return 1; /* failure */
    }
    
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
cmdline_parser_internal (int argc, char * const *argv, struct gengetopt_args_info *args_info,
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
        { "positives",	1, NULL, 'p' },
        { "environment",	1, NULL, 'e' },
        { "data",	1, NULL, 'd' },
        { "genome",	1, NULL, 'g' },
        { "distance",	1, NULL, 'D' },
        { "normalize",	0, NULL, 'N' },
        { "zscore",	0, NULL, 'Z' },
        { "skip",	1, NULL, 'S' },
        { "memmap",	0, NULL, 'm' },
        { "verbosity",	1, NULL, 'v' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "hVp:e:d:g:D:NZS:mv:", long_options, &option_index);

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
        case 'p':	/* Positive gene list.  */
        
        
          if (update_arg( (void *)&(args_info->positives_arg), 
               &(args_info->positives_orig), &(args_info->positives_given),
              &(local_args_info.positives_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "positives", 'p',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* List of environment features and default values.  */
        
        
          if (update_arg( (void *)&(args_info->environment_arg), 
               &(args_info->environment_orig), &(args_info->environment_given),
              &(local_args_info.environment_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "environment", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Feature values for each data set.  */
        
        
          if (update_arg( (void *)&(args_info->data_arg), 
               &(args_info->data_orig), &(args_info->data_given),
              &(local_args_info.data_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "data", 'd',
              additional_error))
            goto failure;
        
          break;
        case 'g':	/* SGD features file.  */
        
        
          if (update_arg( (void *)&(args_info->genome_arg), 
               &(args_info->genome_orig), &(args_info->genome_given),
              &(local_args_info.genome_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "genome", 'g',
              additional_error))
            goto failure;
        
          break;
        case 'D':	/* Similarity measure.  */
        
        
          if (update_arg( (void *)&(args_info->distance_arg), 
               &(args_info->distance_orig), &(args_info->distance_given),
              &(local_args_info.distance_given), optarg, cmdline_parser_distance_values, "pearnorm", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "distance", 'D',
              additional_error))
            goto failure;
        
          break;
        case 'N':	/* Normalize distances.  */
        
        
          if (update_arg((void *)&(args_info->normalize_flag), 0, &(args_info->normalize_given),
              &(local_args_info.normalize_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "normalize", 'N',
              additional_error))
            goto failure;
        
          break;
        case 'Z':	/* Convert correlations to z-scores.  */
        
        
          if (update_arg((void *)&(args_info->zscore_flag), 0, &(args_info->zscore_given),
              &(local_args_info.zscore_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "zscore", 'Z',
              additional_error))
            goto failure;
        
          break;
        case 'S':	/* PCL columns to skip after ID.  */
        
        
          if (update_arg( (void *)&(args_info->skip_arg), 
               &(args_info->skip_orig), &(args_info->skip_given),
              &(local_args_info.skip_given), optarg, 0, "2", ARG_INT,
              check_ambiguity, override, 0, 0,
              "skip", 'S',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* Memory map input DABs.  */
        
        
          if (update_arg((void *)&(args_info->memmap_flag), 0, &(args_info->memmap_given),
              &(local_args_info.memmap_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "memmap", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Message verbosity.  */
        
        
          if (update_arg( (void *)&(args_info->verbosity_arg), 
               &(args_info->verbosity_orig), &(args_info->verbosity_given),
              &(local_args_info.verbosity_given), optarg, 0, "5", ARG_INT,
              check_ambiguity, override, 0, 0,
              "verbosity", 'v',
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
