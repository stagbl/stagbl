#!/usr/bin/env python
#################################################################################
#                       StagBL Configure Script                                 #
#################################################################################
#                                                                               #
# Heavily borrows from HPGMG (https://bitbucket.org/hpgmg/hpgmg)                #
#                                                                               #
# A central concept is that of an "arch", a name for a particular configuration.#
# This defines the name of a subdirectory which (by default) build products     #
# end up in. This allows multiple parallel builds.                              #
#################################################################################
import os
import errno
import sys

try :
    import json
except ImportError:
    print("""ERROR: Could not import json""")
    sys.exit(1)

try:
    import argparse
except ImportError:
    print("""ERROR: Could not import argparse""")
    sys.exit(1)

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise

def main():
    parser = argparse.ArgumentParser(description='Configure StagBL')
    parser.add_argument('--arch', help='Name of this configuration', default=None)
    po = parser.add_argument_group('PETSc')
    po.add_argument('--petsc-dir', help='PETSC_DIR', default=os.environ.get('PETSC_DIR',''))
    po.add_argument('--petsc-arch', help='PETSC_ARCH', default=os.environ.get('PETSC_ARCH',''))
    cf = parser.add_argument_group('Compilers and flags')
    cf.add_argument('--CC', help='Path to C compiler', default=os.environ.get('CC',''))
    cf.add_argument('--CFLAGS', help='Flags for C compiler', default=os.environ.get('CFLAGS',''))
    cf.add_argument('--CPPFLAGS', help='Flags for C preprocessor', default=os.environ.get('CPPFLAGS',''))
    cf.add_argument('--LDFLAGS', help='Flags to pass to linker', default=os.environ.get('LDFLAGS',''))
    cf.add_argument('--LDLIBS', help='Libraries to pass to linker', default=os.environ.get('LDLIBS',''))
    args = parser.parse_args()
    if not args.arch:
        args.arch = args.petsc_arch
    if not args.arch:
        args.arch = 'build'
    mkdir_p(args.arch)
    configure(args)

def configure(args):
    open(os.path.join(args.arch,'variables.mk'), 'w').write(variables(args))
    open(os.path.join(args.arch,'Makefile'), 'w').write(makefile(args))
    dumpinfo(args)
    print('Configuration complete in: %s' % os.path.realpath(args.arch))
    print('To build:')
    print('make -j3 -C %s' % args.arch)
    print('To test:')
    print('make -j3 -C %s test' % args.arch)

def variables(args):
    if args.CC:
        CC = args.CC
    else:
        CC = '$(PCC)'
    m = ['STAGBL_ARCH = %s' % args.arch,
         'STAGBL_CC = %s' % CC,
         'STAGBL_CFLAGS = %s' % (args.CFLAGS if args.CFLAGS else '$(PCC_FLAGS) '),
         'STAGBL_CPPFLAGS = %s' % ('$(CCPPFLAGS) ' + args.CPPFLAGS),
         'STAGBL_LDFLAGS = %s' % args.LDFLAGS,
         'STAGBL_LDLIBS = %s' % args.LDLIBS,
         'SRCDIR = %s' % os.path.abspath(os.path.dirname(__name__)),
         'STAGBL_DIR = %s' % os.path.abspath(os.path.dirname(__name__)), # the same as SRCDIR
         ]
    m.append('STAGBL_LIB = -L$(STAGBL_DIR)/$(STAGBL_ARCH)/lib -lstagbl $(PETSC_LIB)')
    if not args.petsc_dir:
        raise RuntimeError('PETSC_DIR must be provided if using PETSc')
    m.append('PETSC_DIR = %s' % args.petsc_dir)
    m.append('PETSC_ARCH = %s' % args.petsc_arch)
    found = False
    for variables_path in [os.path.join('lib', 'petsc', 'conf', 'variables'),
                           os.path.join('lib', 'petsc-conf', 'variables'),
                           os.path.join('conf', 'variables')]:
        if os.path.exists(os.path.join(args.petsc_dir,variables_path)):
            m.append('include $(PETSC_DIR)/' + variables_path +'\n')
            found = True
    if not found:
        raise RuntimeError('Could not find PETSc variables file in PETSC_DIR=%s' % (args.petsc_dir,))
    return '\n'.join(m)

def makefile(args):
    m = ['include variables.mk']
    m.append('include $(STAGBL_DIR)/stagbl.mk\n')
    return '\n'.join(m)

def dumpinfo(args):
    data = {}
    data['command'] = ' '.join(sys.argv)
    with open(os.path.join(args.arch,'config.json'),'w') as f :
        json.dump(data,f)

if __name__ == "__main__" :
    main()
