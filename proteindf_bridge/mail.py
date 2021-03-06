#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
# 
# This file is a part of the ProteinDF software package.
# 
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import os
import smtplib

from email.mime.text import MIMEText
from email.header import Header
from email.utils import formatdate
try:
    import configparser
except ImportError:
    # for Python 2.x
    import ConfigParser as configparser

class Mail(object):
    def __init__(self):
        self.smtp_server = "localhost"
        self.smtp_port = None
        self.use_SSL = False
        self.smtp_account = ""
        self.smtp_password = ""
        self.from_address = ""
        self.to_address = ""
        self.charset = "ISO-2022-JP"
        self.subject = ""
        self.text = ""

    def load_config(self, path):
        ini = ConfigParser.SafeConfigParser()
        if os.path.exists(path):
            f = open(INI_FILE, "r")
            ini.readfp(f)
            f.close()

        self.smtp_server = ini.get('mail', 'smtp_server')
        self.smtp_port = ini.get('mail', 'smtp_port')
        self.use_SSL = ini.get('mail', 'use_SSL')
        self.smtp_account = ini.get('mail', 'smtp_account')
        self.smtp_password = ini.get('mail', 'smtp_password')
        self.from_address = ini.get('mail', 'from_address')
        
    def save_config(self, path):
        ini = ConfigParser.SafeConfigParser()
        ini.add_section('mail')
        ini.set('mail', 'smtp_server', self.smtp_server)
        ini.set('mail', 'smtp_port', str(self.smtp_port))
        ini.set('mail', 'use_SSL', str(self.use_SSL))
        ini.set('mail', 'smtp_account', self.smtp_account)
        ini.set('mail', 'smtp_password', self.smtp_password)
        ini.set('mail', 'from_address', self.from_address)
        
        f = open(path, 'w')
        ini.write(f)
        f.close()
        
    def send(self):
        msg = MIMEText(self.text.encode(self.charset),
                       "plain",
                       self.charset)
        msg["Subject"] = Header(self.subject, self.charset)
        msg["From"]    = self.from_address
        msg["To"]      = self.to_address
        msg["Date"]    = formatdate(localtime=True)

        smtp = None
        if self.use_SSL:
            if self.smtp_port:
                smtp = smtplib.SMTP_SSL(self.smtp_server,
                                        self.smtp_port)
            else:
                smtp = smtplib.SMTP_SSL(self.smtp_server)
        else:
            if self.smtp_port:
                smtp = smtplib.SMTP(self.smtp_server,
                                    self.smtp_port)
            else:
                smtp = smtplib.SMTP(self.smtp_server)
            smtp.ehlo() # ehlo_or_helo_if_needed()
            smtp.starttls()

        smtp.ehlo()
        smtp.login(self.smtp_account,
                   self.smtp_password)

        smtp.sendmail(self.from_address,
                      self.to_address,
                      msg.as_string())
        smtp.close()

    # smtp_server ==============================================================
    def _get_smtp_server(self):
        return self._smtp_server

    def _set_smtp_server(self, value):
        self._smtp_server = str(value)
        
    smtp_server = property(_get_smtp_server, _set_smtp_server)

    # smtp_port ================================================================
    def _get_smtp_port(self):
        return self._smtp_port

    def _set_smtp_port(self, value):
        self._smtp_port = value
        
    smtp_port = property(_get_smtp_port, _set_smtp_port)

    # use_SSL ==================================================================
    def _get_use_SSL(self):
        return self._use_SSL

    def _set_use_SSL(self, value):
        self._use_SSL = bool(value)

    use_SSL = property(_get_use_SSL, _set_use_SSL)
        
    # smtp_account =============================================================
    def _get_smtp_account(self):
        return self._smtp_account

    def _set_smtp_account(self, value):
        self._smtp_account = str(value)

    smtp_account = property(_get_smtp_account, _set_smtp_account)
    
    # smtp_password ============================================================
    def _get_smtp_password(self):
        return self._smtp_password

    def _set_smtp_password(self, value):
        self._smtp_password = str(value)
    
    smtp_password = property(_get_smtp_password, _set_smtp_password)
    
    # from_address =============================================================
    def _get_from_address(self):
        return self._from_address

    def _set_from_address(self, value):
        self._from_address = str(value)

    from_address = property(_get_from_address, _set_from_address)
        
    # to_address ===============================================================
    def _get_to_address(self):
        return self._to_address

    def _set_to_address(self, value):
        self._to_address = str(value)

    to_address = property(_get_to_address, _set_to_address)

    # charset ==========================================================
    def _get_charset(self):
        return self._charset

    def _set_charset(self, value):
        self._charset = str(value)

    charset = property(_get_charset, _set_charset)
    
    # subject ==========================================================
    def _get_subject(self):
        return self._subject

    def _set_subject(self, subject):
        self._subject = str(subject)

    subject = property(_get_subject, _set_subject)
        
    # text =============================================================
    def _get_text(self):
        return self._text
    
    def _set_text(self, contents):
        self._text = str(contents)

    text = property(_get_text, _set_text)
        
        


if __name__ == "__main__":
    import doctest
    doctest.testmod()
