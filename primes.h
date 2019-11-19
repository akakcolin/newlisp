/* primes.h - table of primitives

    Copyright (C) 2015 Lutz Mueller

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef PRIMES_H
#define PRIMES_H

PRIMITIVE primitive[] =
	{
	/* --------- core ------------------ */
#ifdef SYMBOL_DEBUG
	{"dump-symbol",		p_dumpSymbol,	0},
#endif
	{"while",		p_while,	2},
	{"until",		p_until,	2},
	{"do-while",	p_doWhile,	2},
	{"do-until",	p_doUntil,	2},
	{"dotimes",		p_dotimes,	2},
	{"for",			p_for,		2},
	{"sequence",	p_sequence,	0},
	{"series",		p_series,	0},
	{"doargs",		p_doargs,	2},
	{"dolist",		p_dolist,	2},
	{"dostring",	p_dostring, 2},
	{"dotree",		p_dotree,	2},
	{"case",		p_case,		2},
	{"cond",		p_condition,1},
	{"begin",		p_evalBlock,1}, 
	{"and",			p_and,		0},
	{"if",			p_if,		2},
	{"if-not",		p_ifNot,	2},
	{"unless",		p_unless,	2},
	{"when",		p_when, 	2},
	{"or",			p_or,		0},
	{"quote",		p_quote,	0},
	{"silent",		p_silent,	0},
	{"eval",		p_eval,		0},
	{"amb",			p_amb,		0},
	{"catch",		p_catch,	0},
	{"throw",       p_throw,    0},
	{"apply",		p_apply,	0},
	{"curry",		p_curry,	0},
	{"args",		p_args,		0},
	{"map",			p_map,		0},
	{"term",       	p_term, 	0},
	{"filter",		p_filter,	0},
	{"clean",		p_clean,	0},
	{"index",		p_index,	0},
	{"define",		p_define,	0x402},
	{"define-macro",	p_defineMacro,	0x402},
    {"macro",       p_macro,    2},
	{"expand",		p_expand,	0},
	{"constant",	p_constant,	0x400},
	{"global",		p_global,	0},
	{"set",			p_set,		0x400},
	{"setf",		p_setf,		0x400},
	{"setq",		p_setf,		0x400},
	{"let",			p_let,		0x402},
	{"local",		p_local,	2},
	{"letn",		p_letn,		0x002},
	{"letex",		p_letExpand,0x403},
	{"first",		p_first,	0},
	{"flat",		p_flat,		0},
	{"last",		p_last,		0},
	{"rest",		p_rest,		0},
	{"cons",		p_cons,		0},
	{"append",		p_append,	0},
	{"extend",		p_extend,	0x400},
	{"list",		p_list,		0},
	{"nth",			p_nth,		0},
	{"ref",			p_ref,		0},
	{"ref-all",		p_refAll,	0},
	{"set-ref",		p_setRef,	0x400},
	{"set-ref-all",	p_setRefAll,0x400},
	{"select",		p_select,	0},
	{"collect",		p_collect,	0},
	{"swap",		p_swap,		0}, 
	{"slice",		p_slice,	0},
	{"length",		p_length,	0},
	{"find",		p_find,		0},
	{"search",		p_search,	0},
	{"member",		p_member,	0},
	{"intersect",	p_intersect,	0},
	{"difference",	p_difference,	0},
	{"union",	    p_union,	0},
	{"unique",		p_unique,	0},
	{"assoc",		p_assoc,	0},
	{"lookup",		p_lookup,	0},
	{"count",		p_count,	0},
	{"pop-assoc",	p_popAssoc,	0x400},
	{"replace",		p_replace,	0x400},
	{"sort",		p_sort,		0x400},
	{"push",		p_push,		0x400},
	{"pop",			p_pop,		0x400},
	{"reverse",		p_reverse,	0x400},
	{"rotate",		p_rotate,	0x400},
	{"dup",			p_dup,		0},
	{"not",			p_not,		0},
	{"+",			p_add,		0},
	{"-",			p_subtract,	0},
	{"*",			p_multiply,	0},
	{"/",			p_divide,	0},
	{"%",			p_modulo,	0},
	{"<",			p_less,	0},
	{">",			p_greater,	0},
	{"<=",			p_lessEqual,0},
	{">=",			p_greaterEqual,0},
	{"=",			p_equal,	0},
	{"!=",			p_notEqual,	0},
	{"++",			p_incrementI,0x400},
	{"--",			p_decrementI,0x400},

	/* --------- bit ops --------------- */
	{"<<",			p_shiftLeft,0},
	{">>",			p_shiftRight,	0},
	{"&",			p_bitAnd,	0},
	{"|",			p_bitOr,	0},
	{"^",			p_bitXor,	0},
	{"~",			p_bitNot,	0},

	/* --------- math and float ------- */
	{"inc",		    p_incrementF,0x400},
	{"dec",		    p_decrementF,0x400},
	{"add",		    p_addFloat,	0},
	{"sub",		    p_subFloat,	0},
	{"mul",		    p_mulFloat,	0},
	{"div",		    p_divFloat,	0},
	{"abs",		    p_abs,		0},
	{"ceil",		p_ceil,     0},
	{"floor",		p_floor,	0},
	{"erf",			p_erf,		0},
	{"sgn",			p_sgn,		0},
	{"sin",		    p_sin,		0},
	{"cos",		    p_cos,		0},
	{"tan",		    p_tan,		0},
	{"asin",		p_asin,	    0},
	{"acos",		p_acos,	    0},
	{"atan",		p_atan,	    0},
	{"atan2",		p_atan2,	0},
	{"sinh",		p_sinh,		0},
	{"cosh",		p_cosh,		0},
	{"tanh",		p_tanh,		0},
	{"asinh",		p_asinh,	0},
	{"acosh",		p_acosh,	0},
	{"atanh",		p_atanh,	0},
	{"round",		p_round,	0},
	{"exp",		    p_exp,		0},
	{"log",		    p_log,		0},
	{"sqrt",		p_sqrt,     0},
	{"ssq",		    p_ssq,      0},
	{"rand",		p_rand, 	0},
	{"seed",		p_seed,     0},
	{"random",      p_random,	0},
	{"normal",		p_normal,	0},
	{"randomize",	p_randomize,0},
	{"encrypt",		p_encrypt,	0},
	{"min",		    p_minFloat,	0},
	{"max",		    p_maxFloat,	0},
	{"pow",		    p_powFloat,	0},
	{"mod",		    p_modFloat,	0},
	{"prob-z",		p_probabilityZ, 0},
	{"prob-chi2",	p_probabilityChi2, 0},
	{"prob-t",		p_probabilityT, 0},
	{"prob-f",		p_probabilityF, 0},
	{"crit-chi2",	p_criticalChi2,	0},
	{"crit-z",		p_criticalZ,	0},
	{"crit-t",		p_criticalT,	0},
	{"crit-f",		p_criticalF,	0},
	{"fft",		    p_fft,		0},
	{"ifft",		p_ifft,	        0},
	{"beta",		p_beta,	        0},
	{"betai",		p_betai,	0},
	{"gammaln",		p_gammaln,	0},
	{"gammai",		p_gammai,	0},
	{"binomial",	p_binomial,	0},
	{"factor",	p_factor,	0},
	{"pmt",         p_pmt,          0},
	{"pv",          p_pv,           0},
	{"fv",          p_fv,           0},
	{"nper",        p_nper,         0},
	{"npv",         p_npv,          0},
	{"irr",			p_irr,		0},
	{"transpose",	p_matTranspose, 0},
	{"multiply",	p_matMultiply,  0},
	{"invert",		p_matInvert,    0},
	{"det",			p_determinant,  0},
	{"mat",			p_matScalar,		0},
	{"array",		p_array,	0},
	{"array-list",	p_arrayList,	0},
	{"flt",			p_flt,		0},
	{"bayes-train",	p_bayesTrain,	0},
	{"bayes-query",	p_bayesQuery,	0},
#ifdef KMEANS
    {"kmeans-train",p_kmeansTrain, 0},
    {"kmeans-query",p_kmeansQuery, 0},
#endif
    {"stats",       p_stats,        0},
    {"t-test",      p_ttest,        0},
    {"corr",        p_corr,         0},
	{"unify",		p_unify,		0},
	{"bind",		p_bind,			0x400},
	{"uuid",		p_uuid,			0},
	{"gcd",			p_gcd,			0},
	
	/* ------------ string ops ------------- */
	{"eval-string",	p_evalString,	0},
#ifdef EMSCRIPTEN
    {"eval-string-js", p_evalStringJS, 0},
#endif
	{"read-expr",	p_readExpr,		0},
	{"join",		p_join,	        0},
	{"chop",		p_chop,	        0},
	{"explode",		p_explode,		0},
	{"trim",		p_trim,	        0},
	{"char",		p_char,	        0},
	{"starts-with",	p_startsWith,	0},
	{"ends-with",	p_endsWith,	0},
	{"upper-case",	p_upper,	0},
	{"lower-case",	p_lower,	0},
	{"title-case",	p_title,	0},
	{"format",		p_format,	0},
	{"match",		p_match,	0},
	{"regex",		p_regex,	0},
	{"regex-comp",	p_regexComp,0},
	{"int",			p_integer,	0},
	{"integer",		p_integer,	0},
	{"float",		p_float,	0},
	{"string",		p_string,	0},
	{"bits",		p_bits,		0},
	{"get-float",	p_getFloat,	0},
	{"get-string",	p_getString,	0},
	{"get-int",		p_getInteger,	0},
	{"get-long",	p_getLong,	0},
	{"get-char",	p_getChar,	0},
	{"sym",			p_symbol,	0},
	{"parse",		p_parse,	0},
	{"pack",		p_pack, 	0},
	{"unpack",		p_unpack, 	0},
#ifdef XML_SUPPORT
	{"xml-parse",	p_XMLparse, 	0},
	{"xml-error",	p_XMLerror,	0},
	{"xml-type-tags",   p_XMLtypeTags,  0},
#endif
    {"json-parse",  p_JSONparse, 0},
    {"json-error",  p_JSONerror, 0},
	{"base64-enc",		p_base64Enc,	0},
	{"base64-dec",		p_base64Dec,	0},
	{"crc32",		p_crc32,	0},
	{"find-all",	p_findAll,	0},
#ifdef SUPPORT_UTF8
	{"unicode",		p_unicode,	0},
	{"utf8",		p_utf8,		0},
	{"utf8len",		p_utf8len,	0},
#endif

	/* -------------- I/O ------------------ */
	{"pretty-print",	p_prettyPrint,	0},
	{"print",		p_print,	0},
	{"println",		p_println,	0},
	{"read-line",	p_readLine,	0},
	{"write-line",	p_writeLine,	0},
	{"current-line",p_currentLine,0},
	{"device",		p_device,	0},
	{"load",		p_load,		0},
	{"save",		p_save,		0},
	{"source",		p_symbolSource,	0},
	{"open",		p_open,		0},
	{"seek",		p_seek,		0},
	{"close",		p_close,	0},
	{"read-char",	p_readChar,	0},
#ifdef SUPPORT_UTF8
	{"read-utf8",	p_readUTF8,	0},
#endif
	{"write-char",	p_writeChar,0},
	{"read",		p_readBuffer,0x400},
	{"read-buffer",	p_readBuffer,0x400},
	{"write",		p_writeBuffer,0},
	{"write-buffer",p_writeBuffer,0},
	{"write-file",	p_writeFile,0},
	{"append-file",	p_appendFile,	0},
	{"read-file",	p_readFile, 0},
	{"read-key",	p_readKey,	0},
#ifndef EMSCRIPTEN
	{"get-url",		p_getUrl,	0},
	{"put-url",		p_putUrl, 	0},
	{"post-url",	p_postUrl,	0},
	{"delete-url",	p_deleteUrl,0},
	{"destroy",		p_destroyProcess, 0},
	{"exec",        p_exec,     0},
	{"process",		p_process,	0},
	{"pipe",		p_pipe,		0}, 
#ifndef NO_FORK
	{"fork",		p_fork,		0},
	{"wait-pid",	p_waitpid,	0},
#endif
#ifndef NO_SPAWN
	{"spawn",		p_spawn,	0},
	{"sync",		p_sync,		0},
	{"abort",		p_abort,	0},
	{"send",		p_send,		0},
	{"receive",		p_receive,	0},
#endif
#ifndef NO_SHARE
	{"share",		p_share,	0},
#endif
#ifndef NO_SEMAPHORE
	{"semaphore",	p_semaphore,0},
#endif
#ifndef WINDOWS
	{"peek",		p_peek,		0},
#endif

#endif /* ifndef EMSCRIPTEN */

	
	/* ---------  system --------- */
	{"set-locale",		p_setLocale,	0},
	{"for-all",		p_forAll,	0},
	{"exists",		p_exists,	0},
	{"symbols",		p_symbols,	0},
	{"exit",		p_exit,		0},
#ifdef DEBUGGER
	{"debug",		p_debug,	0},
	{"trace-highlight", 	p_traceHighlight,0},
#endif
	{"trace",		p_trace,	0},
	{"reset",		p_reset,	0},
	{"throw-error",	p_throwError,	0},
	{"error-event",	p_errorEvent,	0},
	{"prompt-event",p_promptEvent,	0},
	{"command-event", p_commandEvent, 0},
	{"xfer-event", 	p_transferEvent, 0},
	{"reader-event", 	p_readerEvent, 0},
	{"last-error", 	p_lastError, 0},
	{"history", 	p_history, 0},

#ifndef EMSCRIPTEN
#ifndef NO_SIGNAL
	{"signal",		p_signal,	0},
#endif
#endif
	{"sys-info",	p_systemInfo,	0},
	{"sys-error",	p_systemError,	0},
	{"!",			p_system,	0},
	{"file-info",	p_fileInfo,	0},
	{"copy",		p_copy,	0},
	{"copy-file",	p_copyFile,	0},
	{"rename-file",	p_renameFile,	0},
	{"delete-file",	p_deleteFile,	0},
	{"make-dir",	p_makeDir,	0},
	{"remove-dir",	p_removeDir,	0},
	{"change-dir",	p_changeDir,	0},
	{"directory",	p_directory,	0},
	{"real-path",	p_realpath,		0},
	{"main-args",	p_mainArgs,	0},
	{"env",			p_env,		0},
	{"context",		p_context,	0},
	{":",			p_colon,	0},
	{"self",		p_self,		0},
	{"prefix",		p_prefix,	0},
	{"default",		p_default,	0},
#ifndef EMSCRIPTEN
	{"timer",       p_timerEvent,	0},
	{"import",		p_importLib,	0},
	{"callback",	p_callback,		0},
#ifdef FFI
    {"struct",      p_struct,       0},
#endif
#endif /* ifndef EMSCRIPTEN */
	{"delete",		p_deleteSymbol,	0},
	{"new",			p_new,		0},
	{"def-new",		p_defineNew,	0},
	{"address",		p_address,	0},
	{"dump",		p_dump,	        0},
	{"cpymem",		p_copyMemory,	0},
	{"sleep",		p_sleep,	0},
	{"$",			p_systemSymbol, 0},

    /* --------- predicates ------ */
	{"nil?",		p_isNil,	0},
	{"true?",		p_isTrue,	0},
	{"NaN?",		p_isnan,	0},
	{"inf?",		p_isinf,	0},
	{"integer?",	p_isInteger,	0},
#ifdef BIGINT
	{"bigint?",	p_isBigInteger,	0},
	{"bigint",	    p_bigInt,	0},
#endif
	{"float?",		p_isFloat,	0},
	{"number?",		p_isNumber,	0},
	{"string?",		p_isString,	0},
	{"symbol?",		p_isSymbol,	0},
	{"legal?",		p_isLegal,	0},
	{"context?",	p_isContext,	0},
	{"primitive?",	p_isPrimitive,	0},
	{"atom?",		p_isAtom,	0},
	{"quote?",		p_isQuote,	0},
	{"list?",		p_isList,	0},
	{"lambda?",		p_isLambda,	0},
	{"macro?",		p_isMacro,	0},
	{"array?",		p_isArray,	0},
	{"empty?",		p_isEmpty,	0},
	{"null?",		p_isNull,	0},
	{"zero?",		p_isZero,	0},
	{"file?",		p_isFile,	0},
	{"directory?",	p_isDirectory,	0},
	{"global?",		p_isGlobal,	0},
	{"protected?",	p_isProtected,	0},
	{"odd?",	    p_isOdd,	0},
	{"even?",	    p_isEven,	0},

	/* ------------ date and time --------- */
	{"date",		p_date,		0},
	{"time",		p_time,	        0},
	{"time-of-day",	p_timeOfDay,    0},
	{"now",			p_now,			0},
#ifndef WINDOWS
	{"date-parse",	p_dateParse,	0},
	{"parse-date",	p_dateParse,	0},
#endif
	{"date-list",	p_dateList,		0},
	{"date-value",  p_dateValue,	0},

	/* ------------ net working ------------ */
#ifndef EMSCRIPTEN
	{"net-close",		p_netClose,	0},
	{"net-service",		p_netService,	0},
	{"net-connect",		p_netConnect,	0},
	{"net-accept",		p_netAccept,	0},
	{"net-local",		p_netLocal,	0},
	{"net-peer",		p_netPeer,	0},
	{"net-ipv",			p_netIpv,	0},
	{"net-lookup",		p_netLookup,	0},
	{"net-receive",		p_netReceive,	0x400},
	{"net-receive-from",p_netReceiveFrom,0},
	{"net-receive-udp",	p_netReceiveUDP,0},
	{"net-send",		p_netSend,	0},
	{"net-send-to",		p_netSendTo,	0},
	{"net-send-udp",	p_netSendUDP,	0},
	{"net-listen",		p_netListen,	0},
#ifndef WINDOWS
	{"net-packet",		p_netPacket,	0},
	{"net-ping",		p_netPing,	0},
#endif
	{"net-peek",		p_netPeek,	0},
	{"net-select",		p_netSelect,	0},
	{"net-sessions",	p_netSessions,	0},
	{"net-eval",		p_netEval,	0},
	{"net-interface",	p_netInterface, 0},
 	{"net-error",		p_netLastError,	0},
#endif

    /* -------------numpy like-----------*/
#ifdef num

    {"arange",                 p_arange, 0},
    {"zeros",                  p_zeros, 0},
    {"eye",                    p_eye, 0},
    {"identity",               p_identity, 0},
    {"ones",                   p_ones, 0},
    {"full",                   p_full, 0},
    {"mgrid",                  p_mgrid, 0}, // np.mgrid

    {"diag",                   p_diag, 0},
    {"tril",                   p_tril, 0}, // lower triangle  of an array
    {"triu",                   p_tril, 0}, // lower triangle  of an array
    {"vander",                 p_vander, 0}, // create vandermonde matrix
    {"mat",                    p_mat, 0},
    {"mat-from-str",           p_bmat, 0},    
    {"double-matrix",           p_double_matrix, 0},
           
    
    // Discrete fourier Transform
    {"hfft",                   p_hfft, 0},
    {"ihfft",                  p_ihfft, 0},
    
    {"fft2",                   p_fft2, 0},
    {"ifft2",                   p_ifft2, 0},

    // linear algebra
    {"dot", p_dot, 0},
    {"inner", p_inner, 0},
    {"outer", p_outer, 0},
    {"vdot", p_vdot, 0},
    {"matmul",                 p_matmul, 0},
    {"kron",                   p_kron, 0},    
    {"matrix_power", p_matrixPower, 0},
    {"tensordot", p_tensordot, 0},
    {"einsum",                 p_einsum, 0},
    {"einsum-path", p_einsumPath, 0},
    {"linalg-svd", p_svd, 0},
    {"linalg-qr",  p_qr, 0},
    {"linalg-cholesky", p_cholesky, 0},
    {"linalg-eig", p_eig, 0},
    {"linalg-eigh", p_eigh, 0},
    {"linalg-eigvals", p_eigvals, 0},
    {"linalg-eigvalsh", p_eigvalsh, 0},
    {"linalg-norm", p_lgNorm, 0},
    {"matrix-rank", p_matRank, 0},
    {"trace", p_trace, 0},
    {"linalg-solve", p_lgSolve, 0},
    {"linalg-lstsq", p_lgLstsq, 0},
    {"linalg-inv", p_lgInv, 0},
    {"linalg-pinv", p_lgPinv, 0},

    
    // mathematical functions
    {"degress", p_degress, 0},
    {"radians", p_radians, 0},
    {"rad2deg", p_rad2deg, 0},
    {"deg2rad", p_deg2rad, 0},
    {"around", p_around, 0},
    {"rounds",  p_round, 0},

    {"lcm",    p_lcm, 0}, // return the lowest common multiple

    // complex numbers
    {"angle", p_angle, 0},
    {"real", p_real, 0},
    {"imag", p_imag, 0},
    {"conj", p_conj, 0},
    {"conjugate", p_conjugate, 0},

    {"cbrt", p_cbrt, 0}, // return cube-root of an array
    {"square", p_square, 0}, // return square of an array
    {"absolute", p_absolute, 0},
    {"fabs", p_fabs, 0},
    // logic functions
    
    {"complex?", p_isComplex, 0},
    {"array_equal", p_arrayEqual, 0},
    {"graycode",   p_graycode, 0}, // return graycode

    
    {"copy-matrix", p_copy_matrix, 0}, // return a fresh copy of tensor m
    {"contiguousp", p_contiguousp, 0}, // return true if the elements of m are contiguous in memory
    
    {"eigen-symm", p_eigen_symm, 0}, // eigenvalues of real symmetric matrix
    {"eigen-symmv", p_eigen_symmv, 0}, // eigenvalues and eigenvectors of real symmetric matrix
    {"eigen-symm", p_eigen_symm, 0}, // eigenvalues of real symmetric matrix
    {"eigen-herm", p_eigen_herm, 0}, // eigenvalues of complex hermitian matrix
    {"eigen-hermv", p_eigen_hermv, 0}, // eigenvalues and eigenvectors of complex hermitian matrix
    
        
#endif


#ifdef PANDAS
    {"read-csv", p_readcsv, 0},
    {"to-csv", p_tocsv, 0},
    {"read-excel", p_readexcel, 0},
    {"to-excel", p_toexcel, 0},
    {"read-msgpack", p_readmsgpack, 0},
    {"to-msgpack", p_tomsgpack, 0},
    
#endif
    
    /* ------------- atom ------------ */
#ifdef ATOM
    {"specie-masss",           p_specieMass, 0},    // atom mass
    {"specie-names",           p_specieName, 0},
    {"specie-numbers",         p_specieNumber, 0},
    {"specie-magmon",          p_specieMagmon, 0},
    {"specie-charge",          p_specieCharge, 0},
    {"specie-position",        p_speciePosition, 0},

    {"chemical-formula",       p_chemicalFormula,0},

    {"unit-vector",            p_unitVector, 0},
    {"get-angle",                   p_getAngle, 0},
    {"structure-distance",     p_structDist, 0},    // get distance two structure Frobenius norm of the spatical distance between all coordiantes
    {"wrap-positions",         p_wrapPositions, 0},
    {"get-atom-layers",        p_getLayers, 0},
    {"find-mic",               p_findMic, 0},  // finde mini-image represetnto f vector(s) D
    {"lattice-vector",         p_latticeVector, 0},
    {"reciprocal-lattice-vector", p_reciprocalCell, 0}
    {"abc-to-lattice",         p_abctolattice, 0},
    {"lattice-to-abc",         p_latticetoabc, 0},
    {"coords",                 p_coords, 0},
    {"pbc",                    p_pbc, 0},

    {"rotate-atoms",           p_rotateAtoms, 0},
    {"get-dihedral",           p_getDihedral, 0},
    {"set-dihedral",           p_setDihedral, 0},
    {"rattle-atoms",           p_rattleAtoms, 0},   // randomly displace atoms
    {"listofformat",           p_listofformat, 0},
    {"listinformat",           p_listinformat, 0},
    {"quicksort",              p_quicksort, 0},
    {"packsort",               p_packsort, 0},

    /* --------------structure---------------*/
    {"create-cell",            p_createCell, 0},
    {"create-struct"          p_createStruc, 0},
    
    /* ------------- space group ----------------- */
    {"sg-num-to-name",         p_sg_numGetName, 0}, // return the name of a space group (Hermann-Mauguin symbol)
    {"sg-num-to-patn",         p_sg_numGetPatn, 0}, // patterson space group number
    {"sg-num-to-symnum",       p_sg_numGetSymnum, 0}, // get number of sym op
    {"sg-num-to-symop",        p_sg_numGetSymop, 0},  // get a specific symm operation string
    {"sg-num-to-symops",       p_sg_numGetSymops, 0}, // get all symmetry operation strings
    {"sg-name-to-num",         p_sg_nameGetNum, 0},
    {"sg-name-to-patn",        p_sg_nameGetNum, 0},
    {"sg-name-to-symnum",      p_sg_nameGetNum, 0},
    {"sg-name-to-symop",       p_sg_nameGetNum, 0},
    {"sg-name-to-symops",      p_sg_nameGetNum, 0},

    /* --------------*/
    /* ------------mode---------------*/
    {"rdf-xyz",               p_rdfxyz, 0 },    // compute rdf
    {"bfgs",                  p_bfgs, 0},
    {"kp-to-monkhorst-pack",  p_kp2MonkhorstPack, 0},
    {"dft-calculator",        p_dftCalculator, 0},
#endif

	{NULL,NULL,0},
};

#endif /* PRIMES_H */
