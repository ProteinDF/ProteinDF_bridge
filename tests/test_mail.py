#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
from proteindf_bridge.mail import Mail


class TestMail(unittest.TestCase):
    def test_save_and_load_config(self):
        mail = Mail()
        mail.smtp_server = "smtp.example.com"
        mail.smtp_port = 587
        mail.use_SSL = True
        mail.smtp_account = "user@example.com"
        mail.smtp_password = "secret_password"
        mail.from_address = "from@example.com"

        with tempfile.NamedTemporaryFile(suffix=".ini", delete=False) as f:
            temp_path = f.name

        try:
            mail.save_config(temp_path)

            loaded_mail = Mail()
            loaded_mail.load_config(temp_path)

            self.assertEqual(loaded_mail.smtp_server, "smtp.example.com")
            self.assertEqual(loaded_mail.smtp_port, 587)
            self.assertTrue(loaded_mail.use_SSL)
            self.assertEqual(loaded_mail.smtp_account, "user@example.com")
            self.assertEqual(loaded_mail.smtp_password, "secret_password")
            self.assertEqual(loaded_mail.from_address, "from@example.com")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


if __name__ == "__main__":
    unittest.main()
