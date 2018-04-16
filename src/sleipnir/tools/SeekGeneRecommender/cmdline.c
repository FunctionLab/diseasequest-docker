/*
  File autogenerated by gengetopt version 2.22
  generated with the following command:
  /memex/qzhu/usr/bin/gengetopt -iSeekGeneRecommender.ggo --default-optional -u -N -e 

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

const char *gengetopt_args_info_purpose = "Preprocess datasets for Seek";

const char *gengetopt_args_info_usage = "Usage: SeekGeneRecommender [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help               Print help and exit",
  "      --version            Print version and exit",
  "\nMode:",
  "  -e, --pcl                PCL mode, suitable for dataset gene variance \n                             calculation  (default=off)",
  "\nPCL mode:",
  "  -V, --pcl_list=filename  PCL list",
  "  -F, --pcl_dir=directory  PCL directory",
  "\nInput:",
  "  -i, --input=filename     Gene mapping file",
  "  -q, --query=filename     Query file",
  "\nOutput:",
  "  -D, --dir_out=directory  Output directory",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
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

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->pcl_given = 0 ;
  args_info->pcl_list_given = 0 ;
  args_info->pcl_dir_given = 0 ;
  args_info->input_given = 0 ;
  args_info->query_given = 0 ;
  args_info->dir_out_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  args_info->pcl_flag = 0;
  args_info->pcl_list_arg = NULL;
  args_info->pcl_list_orig = NULL;
  args_info->pcl_dir_arg = NULL;
  args_info->pcl_dir_orig = NULL;
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->query_arg = NULL;
  args_info->query_orig = NULL;
  args_info->dir_out_arg = NULL;
  args_info->dir_out_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->pcl_help = gengetopt_args_info_help[3] ;
  args_info->pcl_list_help = gengetopt_args_info_help[5] ;
  args_info->pcl_dir_help = gengetopt_args_info_help[6] ;
  args_info->input_help = gengetopt_args_info_help[8] ;
  args_info->query_help = gengetopt_args_info_help[9] ;
  args_info->dir_out_help = gengetopt_args_info_help[11] ;
  
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
  free_string_field (&(args_info->pcl_list_arg));
  free_string_field (&(args_info->pcl_list_orig));
  free_string_field (&(args_info->pcl_dir_arg));
  free_string_field (&(args_info->pcl_dir_orig));
  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->query_arg));
  free_string_field (&(args_info->query_orig));
  free_string_field (&(args_info->dir_out_arg));
  free_string_field (&(args_info->dir_out_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, char *values[])
{
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
  if (args_info->pcl_given)
    write_into_file(outfile, "pcl", 0, 0 );
  if (args_info->pcl_list_given)
    write_into_file(outfile, "pcl_list", args_info->pcl_list_orig, 0);
  if (args_info->pcl_dir_given)
    write_into_file(outfile, "pcl_dir", args_info->pcl_dir_orig, 0);
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->query_given)
    write_into_file(outfile, "query", args_info->query_orig, 0);
  if (args_info->dir_out_given)
    write_into_file(outfile, "dir_out", args_info->dir_out_orig, 0);
  

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
  if (! args_info->input_given)
    {
      fprintf (stderr, "%s: '--input' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->query_given)
    {
      fprintf (stderr, "%s: '--query' ('-q') option required%s\n", prog_name, (additional_error ? additional_error : ""));
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
        { "version",	0, NULL, 0 },
        { "pcl",	0, NULL, 'e' },
        { "pcl_list",	1, NULL, 'V' },
        { "pcl_dir",	1, NULL, 'F' },
        { "input",	1, NULL, 'i' },
        { "query",	1, NULL, 'q' },
        { "dir_out",	1, NULL, 'D' },
        { NULL,	0, NULL, 0 }
      };

      c = getopt_long (argc, argv, "heV:F:i:q:D:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'e':	/* PCL mode, suitable for dataset gene variance calculation.  */
        
        
          if (update_arg((void *)&(args_info->pcl_flag), 0, &(args_info->pcl_given),
              &(local_args_info.pcl_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "pcl", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'V':	/* PCL list.  */
        
        
          if (update_arg( (void *)&(args_info->pcl_list_arg), 
               &(args_info->pcl_list_orig), &(args_info->pcl_list_given),
              &(local_args_info.pcl_list_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "pcl_list", 'V',
              additional_error))
            goto failure;
        
          break;
        case 'F':	/* PCL directory.  */
        
        
          if (update_arg( (void *)&(args_info->pcl_dir_arg), 
               &(args_info->pcl_dir_orig), &(args_info->pcl_dir_given),
              &(local_args_info.pcl_dir_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "pcl_dir", 'F',
              additional_error))
            goto failure;
        
          break;
        case 'i':	/* Gene mapping file.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'q':	/* Query file.  */
        
        
          if (update_arg( (void *)&(args_info->query_arg), 
               &(args_info->query_orig), &(args_info->query_given),
              &(local_args_info.query_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "query", 'q',
              additional_error))
            goto failure;
        
          break;
        case 'D':	/* Output directory.  */
        
        
          if (update_arg( (void *)&(args_info->dir_out_arg), 
               &(args_info->dir_out_orig), &(args_info->dir_out_given),
              &(local_args_info.dir_out_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "dir_out", 'D',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          /* Print version and exit.  */
          if (strcmp (long_options[option_index].name, "version") == 0)
          {
          
          
            if (update_arg( 0 , 
                 0 , &(args_info->version_given),
                &(local_args_info.version_given), optarg, 0, 0, ARG_NO,
                check_ambiguity, override, 0, 0,
                "version", 'V',
                additional_error))
              goto failure;
            cmdline_parser_free (&local_args_info);
            return 0;
          
          }
          
          break;
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
