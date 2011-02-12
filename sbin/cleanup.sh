#!/bin/bash
echo '(cd '$1'; tar --remove-files -cf '$2'.arg.tar *.arg.gz)'|bash
echo '(cd '$1'; tar --remove-files -czf '$2'.mgr.tgz *.mgr)'|bash
echo '(cd '$1'; tar --remove-files -czf '$2'.log.tgz log)'|bash
echo '(cd '$1'; rm -rf *.stats)'|bash
echo '(cd '$1'; rm -rf log)'|bash

