/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name */
#define CMDLINE_PARSER_PACKAGE "Data2Bnt"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * input_arg;	/**< @brief Positive gene list.  */
  char * input_orig;	/**< @brief Positive gene list original value given at command line.  */
  const char *input_help; /**< @brief Positive gene list help description.  */
  char * features_arg;	/**< @brief List of features (nodes) and default values.  */
  char * features_orig;	/**< @brief List of features (nodes) and default values original value given at command line.  */
  const char *features_help; /**< @brief List of features (nodes) and default values help description.  */
  char * data_arg;	/**< @brief Feature values for each data set.  */
  char * data_orig;	/**< @brief Feature values for each data set original value given at command line.  */
  const char *data_help; /**< @brief Feature values for each data set help description.  */
  char * quants_arg;	/**< @brief Quantization file for membership values.  */
  char * quants_orig;	/**< @brief Quantization file for membership values original value given at command line.  */
  const char *quants_help; /**< @brief Quantization file for membership values help description.  */
  char * genome_arg;	/**< @brief SGD features file.  */
  char * genome_orig;	/**< @brief SGD features file original value given at command line.  */
  const char *genome_help; /**< @brief SGD features file help description.  */
  double fraction_arg;	/**< @brief Fraction of genome to cover with default values (default='1').  */
  char * fraction_orig;	/**< @brief Fraction of genome to cover with default values original value given at command line.  */
  const char *fraction_help; /**< @brief Fraction of genome to cover with default values help description.  */
  int sparse_flag;	/**< @brief Output sparse matrix (default=off).  */
  const char *sparse_help; /**< @brief Output sparse matrix help description.  */
  int comments_flag;	/**< @brief Include informational comments (default=off).  */
  const char *comments_help; /**< @brief Include informational comments help description.  */
  int xrff_flag;	/**< @brief Generate XRFF formatted output (default=off).  */
  const char *xrff_help; /**< @brief Generate XRFF formatted output help description.  */
  int weights_flag;	/**< @brief Weight XRFF features (default=off).  */
  const char *weights_help; /**< @brief Weight XRFF features help description.  */
  int verbosity_arg;	/**< @brief Message verbosity (default='5').  */
  char * verbosity_orig;	/**< @brief Message verbosity original value given at command line.  */
  const char *verbosity_help; /**< @brief Message verbosity help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int features_given ;	/**< @brief Whether features was given.  */
  unsigned int data_given ;	/**< @brief Whether data was given.  */
  unsigned int quants_given ;	/**< @brief Whether quants was given.  */
  unsigned int genome_given ;	/**< @brief Whether genome was given.  */
  unsigned int fraction_given ;	/**< @brief Whether fraction was given.  */
  unsigned int sparse_given ;	/**< @brief Whether sparse was given.  */
  unsigned int comments_given ;	/**< @brief Whether comments was given.  */
  unsigned int xrff_given ;	/**< @brief Whether xrff was given.  */
  unsigned int weights_given ;	/**< @brief Whether weights was given.  */
  unsigned int verbosity_given ;	/**< @brief Whether verbosity was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char * const *argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */